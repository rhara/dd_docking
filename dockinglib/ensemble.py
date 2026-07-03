"""Ensemble preparation: turn N raw (co-crystal) receptor PDBs into N
flexible-side-chain-ready PDBQT receptors ("ensemble members"), each with
its own docking box and pocket flexible-residue selection.

`receptor_to_pdbqt` uses the `-f FLEXRES` argument supported by
`mk_prepare_receptor.py`, which emits `<basename>_rigid.pdbqt` +
`<basename>_flex.pdbqt` when flexible residues are given, or plain
`<basename>.pdbqt` when none are.
"""
from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from . import pocket, receptor_prep


@dataclass
class EnsembleMember:
    member_id: str
    receptor_pdb: str
    rigid_pdbqt: str
    flex_pdbqt: Optional[str]
    center: List[float]
    size: List[float]
    flexres: List[str] = field(default_factory=list)


def receptor_to_pdbqt(fixed_pdb: Path, center: Sequence[float], size: Sequence[float],
                      out_basename: Path, *, flexres: str = "",
                      charge_model: str = "gasteiger") -> Tuple[Path, Optional[Path]]:
    """Run Meeko's `mk_prepare_receptor.py` to produce a Vina-ready PDBQT
    receptor (and, if `flexres` is non-empty, a separate rigid/flex pair)
    plus the docking box. Returns (rigid_pdbqt, flex_pdbqt_or_None)."""
    cmd = ["mk_prepare_receptor.py", "--read_pdb", str(fixed_pdb),
           "-o", str(out_basename), "-p", "-v",
           "--box_center", *map(str, center), "--box_size", *map(str, size),
           "--charge_model", charge_model, "-a", "--default_altloc", "A"]
    if flexres:
        cmd += ["-f", flexres]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("mk_prepare_receptor failed:\n" + r.stdout[-1500:] + r.stderr[-1500:])

    if flexres:
        rigid = out_basename.with_name(out_basename.name + "_rigid.pdbqt")
        flex = out_basename.with_name(out_basename.name + "_flex.pdbqt")
        return rigid, flex
    rigid = out_basename.with_name(out_basename.name + ".pdbqt")
    return rigid, None


def prepare_ensemble_member(raw_pdb: Path, member_id: str, out_dir: Path, *,
                            chain: str = "A", lig_resname: str = "",
                            pocket_cutoff: float = 5.0, max_flexres: Optional[int] = 10,
                            box_padding: float = 5.0,
                            charge_model: str = "gasteiger") -> EnsembleMember:
    """Fix one receptor conformation with PDBFixer, derive its docking box
    and pocket flexible residues from its co-crystal ligand, and prepare
    rigid/flex PDBQT files for Vina."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fixed_pdb, ligand_lines, _n_ss = receptor_prep.prepare_receptor_pdb(
        raw_pdb, out_dir / f"{member_id}_fixed.pdb",
        chain=chain, lig_resname=lig_resname,
    )
    if not ligand_lines:
        raise ValueError(
            f"{member_id}: co-crystal ligand ({lig_resname!r}) not found, "
            "cannot determine the docking box/pocket residues"
        )

    center, size = pocket.compute_box(ligand_lines, padding=box_padding)
    ref_coords = pocket.ligand_lines_to_coords(ligand_lines)
    residues = pocket.find_pocket_residues(fixed_pdb, ref_coords, cutoff=pocket_cutoff, max_residues=max_flexres)
    flexres_str = pocket.format_flexres(residues)

    rigid_pdbqt, flex_pdbqt = receptor_to_pdbqt(
        fixed_pdb, center, size, out_dir / member_id,
        flexres=flexres_str, charge_model=charge_model,
    )

    return EnsembleMember(
        member_id=member_id,
        receptor_pdb=str(fixed_pdb),
        rigid_pdbqt=str(rigid_pdbqt),
        flex_pdbqt=str(flex_pdbqt) if flex_pdbqt else None,
        center=list(center),
        size=list(size),
        flexres=[str(r) for r in residues],
    )


def prepare_ensemble(members: Sequence[Tuple[str, Path, str]], out_dir: Path, *,
                     chain: str = "A", pocket_cutoff: float = 5.0, max_flexres: Optional[int] = 10,
                     box_padding: float = 5.0, charge_model: str = "gasteiger") -> List[EnsembleMember]:
    """`members` is a sequence of (member_id, raw_pdb_path, lig_resname) --
    each conformation may have a differently-named co-crystal ligand.
    Prepares each and writes a `manifest.json` under `out_dir`."""
    out_dir = Path(out_dir)
    ensemble = [
        prepare_ensemble_member(
            raw_pdb, member_id, out_dir, chain=chain, lig_resname=lig_resname,
            pocket_cutoff=pocket_cutoff, max_flexres=max_flexres,
            box_padding=box_padding, charge_model=charge_model,
        )
        for member_id, raw_pdb, lig_resname in members
    ]
    save_manifest(out_dir / "manifest.json", ensemble)
    return ensemble


def save_manifest(path: Path, ensemble: Sequence[EnsembleMember]) -> None:
    Path(path).write_text(json.dumps([asdict(m) for m in ensemble], indent=2, ensure_ascii=False))


def load_manifest(path: Path) -> List[EnsembleMember]:
    data = json.loads(Path(path).read_text())
    return [EnsembleMember(**d) for d in data]
