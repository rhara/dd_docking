"""Ligand preparation: read a `.smi` library and embed/prepare each SMILES
as a PDBQT string for Vina.

Adapted from `prepare_ligand_pdbqt()` in `.archives/a2a-vs/vslib/docking.py`;
`.smi` parsing mirrors `omegalib/io_utils.py:read_smi`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass(frozen=True)
class Ligand:
    ligand_id: str
    smiles: str


def read_smi(path: str | Path) -> List[Ligand]:
    """Read a whitespace/tab-separated .smi file: `SMILES <id>` per line."""
    ligands: List[Ligand] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}:{lineno}: expected 'SMILES  ID', got {line!r}")
            smiles, ligand_id = parts[0], parts[1]
            ligands.append(Ligand(ligand_id=ligand_id, smiles=smiles))
    return ligands


def prepare_ligand_pdbqt(smiles: str, seed: int = 0) -> Optional[str]:
    """SMILES -> 3D -> PDBQT string (ETKDGv3 embed + MMFF minimize + Meeko).
    Returns None on embedding/preparation failure."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:  # noqa: BLE001
        pass
    setups = MoleculePreparation()(mol)
    if not setups:
        return None
    pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    return pdbqt if ok else None
