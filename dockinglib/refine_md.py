"""Lightweight induced-fit-style post-processing: short implicit-solvent
OpenMM MD for only the top-ranked docking hits, to see whether a pose (and
the pocket around it) stays put once the receptor is allowed to move.

Uses a CPU-friendly implicit-solvent setup rather than full explicit-solvent
MD, which would be overkill for rescoring a shortlist. Heavy deps (openmm,
openff, mdtraj, meeko) are imported lazily.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import pandas as pd

from .progress import RefineProgress


def _protonate_protein(protein_raw: Path, out_pdb: Path, ph: float = 7.4) -> None:
    """Add hydrogens to a heavy-atom-only fixed receptor PDB (the output of
    `receptor_prep.fix_receptor`), needed before building an MD system."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(protein_raw))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)


def _ligand_from_pose(pose_pdbqt: Path) -> Tuple[Any, Any]:
    """Load a docked pose PDBQT as an OpenFF Molecule (with Gasteiger
    partial charges) plus its RDKit Mol."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from openff.toolkit import Molecule
    from rdkit import Chem

    pmol = PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True)
    rdmol = RDKitMolCreate.from_pdbqt_mol(pmol)[0]
    best = Chem.Mol(rdmol)
    best.RemoveAllConformers()
    best.AddConformer(rdmol.GetConformer(0), assignId=True)
    off = Molecule.from_rdkit(best, allow_undefined_stereo=True)
    try:
        off.assign_partial_charges("gasteiger")
    except Exception:  # noqa: BLE001
        pass
    return off, best


def _build_simulation(protein_pdb: Path, off_lig, implicit: bool,
                      temperature_k: float, timestep_fs: float) -> Tuple[Any, Any]:
    """Build an OpenMM `Simulation` (CPU platform, Langevin middle
    integrator) for the given protein PDB + OpenFF ligand molecule.

    Returns (simulation, modeller). `modeller` holds the combined
    protein+ligand topology/positions used to write out the complex PDB.
    """
    from openmm import LangevinMiddleIntegrator, Platform, app, unit
    from openmm.app import Modeller, PDBFile
    from openmmforcefields.generators import SystemGenerator

    ff = (["amber14-all.xml", "implicit/gbn2.xml"] if implicit else ["amber14-all.xml"])
    ff_kwargs = {"constraints": app.HBonds, "rigidWater": False,
                 "hydrogenMass": 4 * unit.amu}
    sysgen = SystemGenerator(forcefields=ff, small_molecule_forcefield="gaff-2.11",
                             molecules=[off_lig], forcefield_kwargs=ff_kwargs)

    protein = PDBFile(str(protein_pdb))
    modeller = Modeller(protein.topology, protein.positions)
    lig_top = off_lig.to_topology().to_openmm()
    lig_pos = off_lig.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)

    system = sysgen.create_system(modeller.topology)
    integrator = LangevinMiddleIntegrator(temperature_k * unit.kelvin,
                                          1.0 / unit.picosecond,
                                          timestep_fs * unit.femtoseconds)
    platform = Platform.getPlatformByName("CPU")
    sim = app.Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    return sim, modeller


def _summarize_rmsd(rmsd: Any, *, stable_mean_cutoff: float = 3.0,
                    stable_final_cutoff: float = 4.0) -> Dict[str, Any]:
    """Reduce a per-frame ligand heavy-atom RMSD series (Angstrom, as
    returned by `_ligand_rmsd_series`) to summary statistics plus a
    "stable" pose classification.

    A pose is classified as "stable" when its mean RMSD over the whole
    trajectory is below `stable_mean_cutoff` AND its RMSD in the final
    frame is below `stable_final_cutoff` -- i.e. it neither drifted a lot
    on average nor ended up far from its starting pose.

    Pure numeric logic, no OpenMM/mdtraj calls -- split out from
    `refine_one_hit` so it can be unit-tested without a real MD run.
    """
    import numpy as np

    return {
        "n_frames": int(len(rmsd)),
        "rmsd_mean": round(float(np.mean(rmsd)), 2),
        "rmsd_final": round(float(rmsd[-1]), 2),
        "rmsd_max": round(float(np.max(rmsd)), 2),
        "stable": bool(np.mean(rmsd) < stable_mean_cutoff and rmsd[-1] < stable_final_cutoff),
    }


def _ligand_rmsd_series(top_pdb: Path, dcd: Path):
    """Ligand heavy-atom RMSD (frame-0 reference) after protein-CA
    superposition, in Angstrom."""
    import mdtraj as md
    import numpy as np

    traj = md.load(str(dcd), top=str(top_pdb))
    ca = traj.topology.select("protein and name CA")
    lig = traj.topology.select("not protein and not water and element != H")
    if len(lig) == 0:
        lig = traj.topology.select("resname UNK or resname LIG or resname MOL")
    traj.superpose(traj, frame=0, atom_indices=ca)
    ref = traj.xyz[0, lig, :]
    diff = traj.xyz[:, lig, :] - ref
    return np.sqrt((diff ** 2).sum(axis=(1, 2)) / len(lig)) * 10.0  # nm -> Angstrom


def refine_one_hit(name: str, pose_pdbqt: Path, protein_pdb: Path, workdir: Path, *,
                   implicit: bool = True, temperature_k: float = 300.0,
                   timestep_fs: float = 4.0, equil_ps: float = 20.0,
                   prod_ps: float = 100.0, report_ps: float = 2.0) -> Dict[str, Any]:
    """Minimize -> equilibrate -> run a short production MD for one
    protein+pose complex, and summarize ligand pose stability.

    Args:
        name: Identifier for this hit (used for logging and the returned
            record's "name" field; typically the ligand id).
        pose_pdbqt: Path to the docked pose PDBQT (single best pose) to
            load as the ligand.
        protein_pdb: Path to the (already protonated) receptor PDB to
            build the complex with.
        workdir: Directory to write intermediate/output files to
            (`complex.pdb`, `prod.dcd`); created if missing.
        implicit: If True, use GBn2 implicit solvent; if the system build
            fails with implicit solvent, automatically retries once in
            vacuum (the returned record's "implicit" field reflects what
            was actually used).
        temperature_k: Simulation temperature, in kelvin.
        timestep_fs: Integrator timestep, in femtoseconds.
        equil_ps: Equilibration run length, in picoseconds.
        prod_ps: Production run length, in picoseconds.
        report_ps: Trajectory frame-writing interval during production, in
            picoseconds.

    Returns:
        A dict record with keys: "name", "status" ("ok" or "fail"),
        "implicit" (bool, solvent model actually used), "n_frames" (int),
        "rmsd_mean"/"rmsd_final"/"rmsd_max" (ligand heavy-atom RMSD in
        Angstrom), "stable" (bool, rmsd_mean < 3.0 and rmsd_final < 4.0),
        and "e_min_kcal"/"e_final_kcal" (potential energy in kcal/mol
        before minimization and after production, respectively).
    """
    from openmm import app, unit

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rec = {"name": name, "status": "fail", "implicit": implicit}

    off_lig, _ = _ligand_from_pose(pose_pdbqt)

    try:
        sim, modeller = _build_simulation(protein_pdb, off_lig, implicit, temperature_k, timestep_fs)
    except Exception as e:  # noqa: BLE001
        print(f"[MD {name}] implicit solvent system build failed, retrying in vacuum: {e}", flush=True)
        implicit = False
        sim, modeller = _build_simulation(protein_pdb, off_lig, False, temperature_k, timestep_fs)
    rec["implicit"] = implicit

    top_pdb = workdir / "complex.pdb"
    with open(top_pdb, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)

    sim.minimizeEnergy(maxIterations=500)
    e_min = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

    dt_ps = timestep_fs / 1000.0
    sim.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
    sim.step(max(1, int(equil_ps / dt_ps)))

    dcd = workdir / "prod.dcd"
    n_report = max(1, int(report_ps / dt_ps))
    n_prod = int(prod_ps / dt_ps)
    sim.reporters.append(app.DCDReporter(str(dcd), n_report))
    sim.step(max(n_report, n_prod))

    e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

    rmsd = _ligand_rmsd_series(top_pdb, dcd)
    rec.update({"status": "ok", **_summarize_rmsd(rmsd)})
    rec.update({
        "e_min_kcal": round(float(e_min), 1),
        "e_final_kcal": round(float(e_final), 1),
    })
    return rec


def refine_top_hits(ranked_csv: str, out_dir: str, *, top_n: int = 5,
                    implicit: bool = True, temperature_k: float = 300.0,
                    timestep_fs: float = 4.0, equil_ps: float = 20.0,
                    prod_ps: float = 100.0, report_ps: float = 2.0,
                    show_progress: bool = True) -> pd.DataFrame:
    """Read a `dockinglib-dock` ranked CSV (must have `pose_pdbqt` and
    `receptor_pdb` columns) and MD-refine/rescore its top `top_n` hits.
    Writes and returns `md_rescore.csv`: Vina affinity re-ranked by
    (stable first, then affinity ascending).

    Args:
        ranked_csv: Path to the ranked results CSV produced by
            `screening.write_results_csv` (needs `ligand_id`,
            `best_affinity`, `pose_pdbqt`, and `receptor_pdb` columns).
        out_dir: Output directory for per-hit MD workdirs, protonated
            receptor PDBs, and the final `md_rescore.csv`.
        top_n: Number of top-ranked rows (by input CSV order) to refine.
        implicit: If True, use GBn2 implicit solvent for each hit (falls
            back to vacuum per-hit if the implicit system build fails);
            if False, always run in vacuum.
        temperature_k: Simulation temperature, in kelvin.
        timestep_fs: Integrator timestep, in femtoseconds.
        equil_ps: Equilibration run length per hit, in picoseconds.
        prod_ps: Production run length per hit, in picoseconds.
        report_ps: Trajectory frame-writing interval during production,
            in picoseconds.
        show_progress: If True, print a progress line per completed hit.

    Returns:
        The results DataFrame (also written to `out_dir/md_rescore.csv`),
        with one row per hit including `best_affinity` and, when the MD
        run succeeded, the stability/RMSD/energy fields documented in
        `refine_one_hit`.
    """
    import re

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(ranked_csv)
    df = df[df["pose_pdbqt"].notna() & df["receptor_pdb"].notna()].head(top_n)
    if df.empty:
        raise ValueError(f"{ranked_csv}: pose_pdbqt/receptor_pdb columns are empty")

    protonated_by_receptor: Dict[str, Path] = {}
    progress = RefineProgress(enabled=show_progress)
    safe_re = re.compile(r"[^A-Za-z0-9._-]+")
    results = []
    for rank_i, (_, r) in enumerate(df.iterrows(), start=1):
        receptor_raw = Path(r["receptor_pdb"])
        if receptor_raw not in protonated_by_receptor:
            protonated = out_dir / f"{receptor_raw.stem}_H.pdb"
            if not protonated.exists():
                _protonate_protein(receptor_raw, protonated)
            protonated_by_receptor[receptor_raw] = protonated
        protein_pdb = protonated_by_receptor[receptor_raw]

        safe = f"{rank_i:02d}_" + safe_re.sub("_", str(r["ligand_id"]))[:40]
        workdir = out_dir / safe
        try:
            res = refine_one_hit(
                str(r["ligand_id"]), Path(r["pose_pdbqt"]), protein_pdb, workdir,
                implicit=implicit, temperature_k=temperature_k, timestep_fs=timestep_fs,
                equil_ps=equil_ps, prod_ps=prod_ps, report_ps=report_ps,
            )
        except Exception as e:  # noqa: BLE001
            res = {"name": str(r["ligand_id"]), "status": f"error:{type(e).__name__}"}
        res["best_affinity"] = float(r["best_affinity"])
        results.append(res)
        if res.get("status") == "ok":
            progress.update(res["name"], res["implicit"], res["rmsd_mean"], res["rmsd_final"], res["stable"])

    out = pd.DataFrame(results)
    if "stable" in out:
        out["stable"] = out["stable"].fillna(False)
        out = out.sort_values(["stable", "best_affinity"], ascending=[False, True])
    out_csv = out_dir / "md_rescore.csv"
    out.to_csv(out_csv, index=False)
    return out
