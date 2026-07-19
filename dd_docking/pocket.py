"""Pocket definition: docking box geometry and flexible-residue selection.

`compute_box` computes docking-box center/size from ligand coordinates;
residue-within-cutoff detection and Meeko flexres formatting live here too.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

Coord = Tuple[float, float, float]


def compute_box(ligand_lines: Sequence[str], padding: float = 5.0, *,
                extra_coords: Optional[Sequence[Coord]] = None):
    """Docking box center/size (each axis) spanning the given ligand
    ATOM/HETATM lines, plus `padding` Angstrom on every side.

    `extra_coords` (e.g. flexible side-chain atom coordinates from
    `residue_coords`, see `ensemble.py`) are folded into the same bounding
    box before padding, so a flexible residue reaching further from the
    ligand than a uniform padding would cover doesn't get its movable atoms
    pushed outside the search box (Vina then reports "no conformations
    completely within the search space" and returns no poses at all)."""
    xs, ys, zs = [], [], []
    for ln in ligand_lines:
        if ln[:6] in ("ATOM  ", "HETATM"):
            xs.append(float(ln[30:38]))
            ys.append(float(ln[38:46]))
            zs.append(float(ln[46:54]))
    if not xs:
        raise ValueError("compute_box: no ligand atoms found")
    for x, y, z in (extra_coords or ()):
        xs.append(x)
        ys.append(y)
        zs.append(z)
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
    min_d2: Dict[Residue, float] = {}
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


def residue_coords(pdb_text: str, residues: Iterable[Residue]) -> List[Coord]:
    """Atom coordinates belonging to the given residues (e.g. chosen
    flexible side chains) -- see `compute_box`'s `extra_coords`."""
    wanted = {(r.chain, r.resnum) for r in residues}
    return [coord for chain, resnum, coord in _parse_atom_lines(pdb_text)
            if (chain, resnum) in wanted]


def format_flexres(residues: Sequence[Residue]) -> str:
    """Format residues for Meeko's `mk_prepare_receptor.py -f FLEXRES`,
    e.g. `A:42,A:87,B:23` (chain prefix always given explicitly per
    residue -- unambiguous regardless of chain interleaving)."""
    return ",".join(str(r) for r in residues)
