"""Lightweight induced-fit-style post-processing: short implicit-solvent
OpenMM MD for only the top-ranked docking hits, to see whether a pose (and
the pocket around it) stays put once the receptor is allowed to move.

Adapted from the CPU-friendly implicit-solvent variant in
`.archives/mpro/.archives/vspipe/md.py` (`_protonate_protein`,
`_build_simulation`, `_run_one_md`, `_ligand_rmsd_series`) -- not the full
explicit-solvent `step7_md.py`, which is overkill for rescoring a shortlist.
Heavy deps (openmm, openff, mdtraj, meeko) are imported lazily.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

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


def _ligand_from_pose(pose_pdbqt: Path):
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
                      temperature_k: float, timestep_fs: float):
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
                   prod_ps: float = 100.0, report_ps: float = 2.0) -> dict:
    """Minimize -> equilibrate -> run a short production MD for one
    protein+pose complex, and summarize ligand pose stability."""
    import numpy as np
    from openmm import app, unit

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rec = {"name": name, "status": "fail", "implicit": implicit}

    off_lig, _ = _ligand_from_pose(pose_pdbqt)

    try:
        sim, modeller = _build_simulation(protein_pdb, off_lig, implicit, temperature_k, timestep_fs)
    except Exception as e:  # noqa: BLE001
        print(f"[MD {name}] implicit solvent 系の構築失敗、真空で再試行: {e}", flush=True)
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
    rec.update({
        "status": "ok",
        "n_frames": int(len(rmsd)),
        "rmsd_mean": round(float(np.mean(rmsd)), 2),
        "rmsd_final": round(float(rmsd[-1]), 2),
        "rmsd_max": round(float(np.max(rmsd)), 2),
        "stable": bool(np.mean(rmsd) < 3.0 and rmsd[-1] < 4.0),
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
    """
    import re

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(ranked_csv)
    df = df[df["pose_pdbqt"].notna() & df["receptor_pdb"].notna()].head(top_n)
    if df.empty:
        raise ValueError(f"{ranked_csv}: pose_pdbqt/receptor_pdb 列が空です")

    protonated_by_receptor: dict[str, Path] = {}
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
