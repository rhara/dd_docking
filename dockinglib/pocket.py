"""Pocket definition: docking box geometry and flexible-residue selection.

No existing project has equivalent code (confirmed during exploration) --
`compute_box` is ported from `.archives/a2a-vs/vslib/docking.py`, and
residue-within-cutoff detection / Meeko flexres formatting is new.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

Coord = Tuple[float, float, float]


def compute_box(ligand_lines: Sequence[str], padding: float = 5.0):
    """Docking box center/size (each axis) spanning the given ligand
    ATOM/HETATM lines plus `padding` Angstrom on every side."""
    xs, ys, zs = [], [], []
    for ln in ligand_lines:
        if ln[:6] in ("ATOM  ", "HETATM"):
            xs.append(float(ln[30:38]))
            ys.append(float(ln[38:46]))
            zs.append(float(ln[46:54]))
    if not xs:
        raise ValueError("compute_box: リガンド原子が見つかりません")
    center = [round((min(v) + max(v)) / 2, 3) for v in (xs, ys, zs)]
    size = [round((max(v) - min(v)) + 2 * padding, 3) for v in (xs, ys, zs)]
    return center, size


def _parse_atom_lines(pdb_text: str):
    """Yield (chain, resnum, coord) for every ATOM/HETATM line."""
    for ln in pdb_text.splitlines():
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        chain = ln[21]
        try:
            resnum = int(ln[22:26])
            coord = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        yield chain, resnum, coord


def ligand_lines_to_coords(ligand_lines: Sequence[str]) -> List[Coord]:
    return [c for _, _, c in _parse_atom_lines("\n".join(ligand_lines))]


@dataclass(frozen=True)
class Residue:
    chain: str
    resnum: int

    def __str__(self) -> str:
        return f"{self.chain}:{self.resnum}"


def find_pocket_residues(receptor_pdb: Path, ref_coords: Iterable[Coord],
                         cutoff: float = 5.0, max_residues: Optional[int] = 10) -> List[Residue]:
    """Residues in `receptor_pdb` with at least one atom within `cutoff`
    Angstrom of any of `ref_coords` (typically the co-crystal/reference
    ligand's atom coordinates), kept as Vina flexible side chains.

    `max_residues` caps how many are returned (closest-to-ligand first,
    then re-sorted by (chain, resnum) for readability) -- unlike rigid
    docking, each flexible residue adds real rotatable-bond degrees of
    freedom to Vina's search rather than a precomputed grid, so a whole
    contact shell (easily 20+ residues at a 5 A cutoff on a normal-sized
    ligand) makes docking impractically slow. Vina's own guidance is to
    keep flexible side chains to a small handful; pass `max_residues=None`
    to disable the cap if you really want every contact residue flexible.
    """
    ref = list(ref_coords)
    if not ref:
        return []
    cutoff2 = cutoff * cutoff
    min_d2: dict[Residue, float] = {}
    for chain, resnum, coord in _parse_atom_lines(Path(receptor_pdb).read_text()):
        best = None
        for rx, ry, rz in ref:
            d2 = (coord[0] - rx) ** 2 + (coord[1] - ry) ** 2 + (coord[2] - rz) ** 2
            if d2 <= cutoff2 and (best is None or d2 < best):
                best = d2
        if best is not None:
            residue = Residue(chain, resnum)
            if residue not in min_d2 or best < min_d2[residue]:
                min_d2[residue] = best

    closest = sorted(min_d2, key=lambda r: min_d2[r])
    if max_residues is not None:
        closest = closest[:max_residues]
    return sorted(closest, key=lambda r: (r.chain, r.resnum))


def format_flexres(residues: Sequence[Residue]) -> str:
    """Format residues for Meeko's `mk_prepare_receptor.py -f FLEXRES`,
    e.g. `A:42,A:87,B:23` (chain prefix always given explicitly per
    residue -- unambiguous regardless of chain interleaving)."""
    return ",".join(str(r) for r in residues)
