"""CPU docking via QuickVina2 (`qvina2`, the conda-forge `qvina` package),
a CPU-tuned fork of AutoDock Vina 1.1.2. Supports flexible-side-chain
receptors via `--flex`, same as upstream Vina.

Runs `qvina2` as a subprocess per ligand rather than through a persistent
grid-map object -- there's nothing to reuse across ligands anyway (each
call recomputes the grid map internally), and this mirrors
`gpu_backend.dock_ligand_gpu`'s shape so `screening.py` can treat both
backends the same way.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Sequence, Tuple

from . import io

BIN_NAME = "qvina2"


def dock_ligand(
    rigid_pdbqt: str,
    ligand_pdbqt: str,
    center: Sequence[float],
    size: Sequence[float],
    *,
    flex_pdbqt: Optional[str] = None,
    seed: int = 0,
    exhaustiveness: int = 16,
    n_poses: int = 5,
    cpu: int = 0,
) -> Optional[Tuple[float, str]]:
    """Dock one ligand with QuickVina2. Returns (best_affinity_kcal_mol,
    poses_pdbqt) or None on failure. `poses_pdbqt` includes flexible
    side-chain coordinates when the receptor was set up with a flex PDBQT.
    `cpu` <= 0 (default) lets qvina2 auto-detect all available cores, same
    convention as the old `vina` Python package's `cpu=0`.
    """
    binary = shutil.which(BIN_NAME)
    if binary is None:
        raise RuntimeError(
            f"{BIN_NAME} not found on PATH -- install with: "
            "mamba install -n dd_docking -c conda-forge qvina"
        )

    with tempfile.TemporaryDirectory(prefix="dd_docking_cpu_") as tmp:
        tmp_path = Path(tmp)
        ligand_path = tmp_path / "ligand.pdbqt"
        ligand_path.write_text(ligand_pdbqt)
        out_path = tmp_path / "out.pdbqt"

        args = [
            binary,
            "--receptor", str(rigid_pdbqt),
            "--ligand", str(ligand_path),
            "--out", str(out_path),
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses),
            "--seed", str(seed),
        ]
        if flex_pdbqt:
            args += ["--flex", str(flex_pdbqt)]
        if cpu > 0:
            args += ["--cpu", str(cpu)]

        try:
            subprocess.run(args, capture_output=True, text=True, check=True)
        except (subprocess.CalledProcessError, OSError):
            return None

        if not out_path.exists():
            return None
        poses_pdbqt = out_path.read_text()

    affinity = io.parse_vina_affinity(poses_pdbqt)
    if affinity is None:
        return None
    return affinity, poses_pdbqt
