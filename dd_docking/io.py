"""PDBQT<->RDKit conversion, pose SDF writing, and results CSV I/O.

Pose SDF properties use `affinity` for the docking score (matches one of
`dd_viewer.scoring.SCORE_PROPERTY_CANDIDATES`, so results written here
load directly into the companion protein-ligand viewer, `dd_viewer`).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd


def parse_vina_affinity(pdbqt: str) -> Optional[float]:
    """Extract the best-mode affinity (kcal/mol) from a Vina-family (Vina,
    QuickVina2, Vina-GPU+...) pose PDBQT's "REMARK VINA RESULT:" line.
    Returns None if not found.
    """
    for line in pdbqt.splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split(":", 1)[1].split()[0])
    return None


def pdbqt_string_to_mol(pdbqt: str):
    """Convert a Vina pose PDBQT string (first/best pose) to an RDKit Mol,
    including flexible side-chain atoms when present. Returns None on
    failure.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    try:
        pmol = PDBQTMolecule(pdbqt, skip_typing=True)
        rdmol = RDKitMolCreate.from_pdbqt_mol(pmol)[0]
    except Exception:  # noqa: BLE001
        return None
    best = Chem.Mol(rdmol)
    best.RemoveAllConformers()
    best.AddConformer(rdmol.GetConformer(0), assignId=True)
    return best


def write_pose_sdf(out_sdf: str, records: Sequence[Tuple[Any, Dict[str, Any]]]) -> None:
    """Write an SDF where each record is (rdkit_mol, properties_dict)."""
    from rdkit import Chem

    writer = Chem.SDWriter(str(out_sdf))
    for mol, props in records:
        for key, value in props.items():
            mol.SetProp(key, str(value))
        writer.write(mol)
    writer.close()


def write_results_csv(out_csv: str, rows: List[Dict[str, Any]]) -> None:
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def read_results_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path)
