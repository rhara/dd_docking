"""Ensemble virtual screening: dock each ligand against every ensemble
member (receptor conformation) and rank ligands by their best-of-ensemble
affinity.

Shape mirrors `overlaylib/screening.py`'s "run N items, keep the best, rank"
pattern, but here "N" is ensemble members per ligand rather than conformers
per molecule. Parallelizes over the flat (ligand, member) task grid (each
task independently docks one ligand into one member's grid maps), matching
the per-task independence assumption of `parallel_map`.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from . import docking, io, ligand_prep
from .ensemble import EnsembleMember
from .ligand_prep import Ligand
from .parallel import parallel_map
from .progress import DockProgress


@dataclass
class EnsembleHit:
    ligand_id: str
    smiles: str
    best_member: Optional[str]
    best_affinity: Optional[float]
    member_affinities: Dict[str, Optional[float]] = field(default_factory=dict)
    best_pose_pdbqt: Optional[str] = None


def rank_key(hit: EnsembleHit):
    """Sort key for best-of-ensemble ranking: lowest affinity first, with
    ligands that failed on every member (affinity=None) sorted last."""
    return (hit.best_affinity is None, hit.best_affinity)


def _dock_task(task) -> Tuple[str, str, Optional[float], Optional[str]]:
    ligand_id, smiles, member, exhaustiveness, n_poses, seed, cpu = task
    pdbqt = ligand_prep.prepare_ligand_pdbqt(smiles, seed=seed)
    if pdbqt is None:
        return ligand_id, member.member_id, None, None
    v = docking.make_vina(
        member.rigid_pdbqt, member.center, member.size,
        flex_pdbqt=member.flex_pdbqt, seed=seed, exhaustiveness=exhaustiveness, cpu=cpu,
    )
    result = docking.dock_ligand(v, pdbqt, n_poses=n_poses)
    if result is None:
        return ligand_id, member.member_id, None, None
    affinity, poses_pdbqt = result
    return ligand_id, member.member_id, affinity, poses_pdbqt


def screen_ensemble(
    ensemble: Sequence[EnsembleMember],
    ligands: Sequence[Ligand],
    *,
    exhaustiveness: int = 16,
    n_poses: int = 5,
    seed: int = 0,
    n_jobs: int = 1,
    show_progress: bool = True,
    out_csv: Optional[str] = None,
    out_sdf: Optional[str] = None,
    top_n: Optional[int] = None,
) -> List[EnsembleHit]:
    """Dock every ligand against every ensemble member; rank ligands by the
    best (lowest) affinity achieved across the ensemble.

    `n_jobs` parallelizes across the full (ligand, member) task grid: 1
    (default) runs sequentially, letting Vina use every core internally for
    each task; <=0 or N>1 runs that many worker *processes*, each pinned to
    (available cores / worker count) so the workers don't oversubscribe the
    machine and fight each other for cores. Flexible-side-chain docking is
    heavy enough per task (each flexible residue adds real search
    dimensions -- see `pocket.find_pocket_residues`) that Vina's own
    internal multi-threading matters a lot; unlike cheap rigid-receptor
    docking pipelines (which typically pin every job to 1 CPU and rely
    purely on process-level parallelism), pinning every worker to 1 CPU
    here was measured to be ~5x slower per task than letting it use
    several cores.
    """
    if not ensemble:
        raise ValueError("screen_ensemble: ensemble is empty")

    if n_jobs == 1:
        cpu_per_task = 0
    else:
        n_workers = n_jobs if n_jobs and n_jobs > 0 else (os.cpu_count() or 1)
        cpu_per_task = max(1, (os.cpu_count() or 1) // n_workers)
    tasks = [
        (ligand.ligand_id, ligand.smiles, member, exhaustiveness, n_poses, seed, cpu_per_task)
        for ligand in ligands
        for member in ensemble
    ]

    pending: Dict[str, EnsembleHit] = {
        ligand.ligand_id: EnsembleHit(ligand.ligand_id, ligand.smiles, None, None)
        for ligand in ligands
    }
    remaining = {ligand.ligand_id: len(ensemble) for ligand in ligands}

    progress = DockProgress(enabled=show_progress)
    hits: List[EnsembleHit] = []
    for ligand_id, member_id, affinity, poses_pdbqt in parallel_map(
        _dock_task, tasks, n_jobs=n_jobs, announce=show_progress,
    ):
        hit = pending[ligand_id]
        hit.member_affinities[member_id] = affinity
        if affinity is not None and (hit.best_affinity is None or affinity < hit.best_affinity):
            hit.best_affinity = affinity
            hit.best_member = member_id
            hit.best_pose_pdbqt = poses_pdbqt

        remaining[ligand_id] -= 1
        if remaining[ligand_id] == 0:
            hits.append(hit)
            if hit.best_affinity is not None:
                progress.update(hit.ligand_id, hit.best_member, hit.best_affinity)

    hits.sort(key=rank_key)
    if top_n is not None:
        hits = hits[:top_n]

    if out_csv:
        write_results_csv(out_csv, hits, ensemble=ensemble, pose_dir=_pose_dir_for(out_csv))
    if out_sdf:
        write_hits_sdf(out_sdf, hits)

    return hits


def _pose_dir_for(out_csv: str):
    from pathlib import Path
    return Path(out_csv).with_suffix("").parent / (Path(out_csv).stem + "_poses")


def write_results_csv(out_csv: str, hits: Sequence[EnsembleHit], *,
                      ensemble: Optional[Sequence[EnsembleMember]] = None,
                      pose_dir=None) -> None:
    """Write the ranked results table. When `ensemble` is given, each row
    also gets the winning member's fixed receptor PDB path (`receptor_pdb`)
    so `refine_md.py` can rebuild the complex later. When `pose_dir` is
    given, each hit's best pose PDBQT is saved to its own file there and
    referenced by the `pose_pdbqt` column (also required by `refine_md.py`).
    """
    from pathlib import Path

    receptor_by_member = {m.member_id: m.receptor_pdb for m in ensemble} if ensemble else {}
    if pose_dir:
        Path(pose_dir).mkdir(parents=True, exist_ok=True)

    rows = []
    for rank, hit in enumerate(hits, start=1):
        pose_path = None
        if pose_dir and hit.best_pose_pdbqt:
            pose_path = Path(pose_dir) / f"{rank:03d}_{hit.ligand_id}.pdbqt"
            pose_path.write_text(hit.best_pose_pdbqt)
        rows.append({
            "rank": rank,
            "ligand_id": hit.ligand_id,
            "smiles": hit.smiles,
            "best_member": hit.best_member,
            "best_affinity": hit.best_affinity,
            "receptor_pdb": receptor_by_member.get(hit.best_member) if hit.best_member else None,
            "pose_pdbqt": str(pose_path) if pose_path else None,
            **{f"affinity[{m}]": a for m, a in hit.member_affinities.items()},
        })
    io.write_results_csv(out_csv, rows)


def write_hits_sdf(out_sdf: str, hits: Sequence[EnsembleHit]) -> None:
    """Write the best pose of each hit, tagged for the `plviewerlib` viewer
    (`affinity` property is one of its recognized score keys)."""
    records = []
    for rank, hit in enumerate(hits, start=1):
        if not hit.best_pose_pdbqt:
            continue
        mol = io.pdbqt_string_to_mol(hit.best_pose_pdbqt)
        if mol is None:
            continue
        props = {
            "Role": "hit",
            "Rank": rank,
            "LigandId": hit.ligand_id,
            "BestMember": hit.best_member,
            "affinity": f"{hit.best_affinity:.3f}",
        }
        records.append((mol, props))
    io.write_pose_sdf(out_sdf, records)
