"""Thin AutoDock Vina wrapper, extended for flexible-side-chain receptors.

Uses the `flex_pdbqt` argument that Vina's Python API (1.2.7+) supports
(`Vina.set_receptor(rigid_pdbqt_filename, flex_pdbqt_filename)`).
"""
from __future__ import annotations

from typing import Optional, Sequence


def make_vina(rigid_pdbqt: str, center: Sequence[float], size: Sequence[float], *,
             flex_pdbqt: Optional[str] = None, seed: int = 0,
             exhaustiveness: int = 16, cpu: int = 0, verbosity: int = 0):
    """Build a Vina instance with the receptor (and optional flexible side
    chains) and grid maps set up, ready to dock ligands against
    repeatedly. `exhaustiveness` is stashed on the instance for `dock_ligand`.
    """
    from vina import Vina

    v = Vina(sf_name="vina", seed=seed, verbosity=verbosity, cpu=cpu)
    if flex_pdbqt:
        v.set_receptor(rigid_pdbqt, flex_pdbqt)
    else:
        v.set_receptor(rigid_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v._exhaustiveness = exhaustiveness
    return v


def dock_ligand(v, ligand_pdbqt: str, n_poses: int = 5):
    """Dock one ligand PDBQT string against an already-configured Vina
    instance. Returns (best_affinity_kcal_mol, poses_pdbqt) or None on
    failure. `poses_pdbqt` includes flexible side-chain coordinates when
    the receptor was set up with a flex PDBQT.
    """
    try:
        v.set_ligand_from_string(ligand_pdbqt)
        v.dock(exhaustiveness=getattr(v, "_exhaustiveness", 16), n_poses=n_poses)
        return float(v.energies(n_poses=1)[0][0]), v.poses(n_poses=n_poses)
    except Exception:  # noqa: BLE001
        return None
