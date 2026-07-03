"""Receptor structure preparation: fetch/split a raw PDB, tidy it for
downstream tools, and fix missing atoms with PDBFixer.

Heavy dependencies (openmm, pdbfixer) are imported lazily inside functions
so this module stays importable without those installed.

Hydrogens are deliberately NOT added here: Meeko's `mk_prepare_receptor.py`
does its own protonation/charge assignment during PDBQT conversion (see
`ensemble.py`), and pre-adding hydrogens with PDBFixer has been observed to
clash with that step.
"""
from __future__ import annotations

import math
import urllib.request
from pathlib import Path

RCSB_PDB = "https://files.rcsb.org/download/{pdb}.pdb"


def download_pdb(pdb_id: str, dest: Path) -> str:
    """Fetch a raw PDB entry from RCSB and return its text contents.

    Convenience utility for users who want a co-crystal structure straight
    from RCSB before running it through `prepare_receptor_pdb`. Not used
    internally elsewhere in this package (`prepare_receptor_pdb` and
    friends take an already-downloaded PDB path), but kept as public API
    for callers doing their own PDB fetching.

    Args:
        pdb_id: 4-character RCSB PDB accession code (e.g. "6W63").
        dest: Local file path to save the downloaded PDB to. If this path
            already exists, the download is skipped and the cached file's
            text is returned instead (simple on-disk cache).

    Returns:
        The full text contents of the (possibly cached) PDB file.
    """
    dest = Path(dest)
    if not dest.exists():
        urllib.request.urlretrieve(RCSB_PDB.format(pdb=pdb_id), dest)
    return dest.read_text()


def _altloc_ok(line: str) -> bool:
    a = line[16] if len(line) > 16 else " "
    return a in (" ", "A")


def split_receptor(text: str, chain: str, lig_resname: str,
                   drop_resseq_above: int | None = None):
    """Split raw PDB text into (receptor ATOM/TER lines, co-crystal ligand
    HETATM lines) for one chain, keeping only the primary altloc."""
    receptor, ligand, lig_chain = [], [], None
    for ln in text.splitlines():
        rec = ln[:6]
        if rec == "ATOM  " and ln[21] == chain and _altloc_ok(ln):
            if drop_resseq_above is not None:
                try:
                    if int(ln[22:26]) > drop_resseq_above:
                        continue
                except ValueError:
                    pass
            receptor.append(ln)
        elif rec == "HETATM":
            resn = ln[17:20].strip().upper()
            if resn == lig_resname.upper() and _altloc_ok(ln):
                if lig_chain is None:
                    lig_chain = ln[21]
                if ln[21] == lig_chain:
                    ligand.append(ln)
        elif rec == "TER   " and ln[21:22] == chain:
            receptor.append(ln)
    return receptor, ligand


def tidy_receptor(receptor_lines, out_pdb: Path, ss_cut: float = 2.5) -> int:
    """Renumber residues to 1..N, insert TER at backbone breaks, and rename
    disulfide-bonded CYS to CYX (auto-detected by SG-SG distance) so that
    PDBFixer/Meeko don't choke on numbering gaps or miss disulfides.
    Returns the number of disulfide bonds found.
    """
    lines = list(receptor_lines)
    new_idx, cur, prev = {}, 0, None
    sg, bbN, bbC = [], {}, {}
    for n, ln in enumerate(lines):
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        rid = (ln[21], ln[22:26], ln[26])
        if rid != prev:
            cur += 1
            prev = rid
        new_idx[n] = cur
        name = ln[12:16].strip()
        crd = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if name == "N":
            bbN[cur] = crd
        elif name == "C":
            bbC[cur] = crd
        elif name == "SG" and ln[17:20].strip() == "CYS":
            sg.append((cur, crd))
    ss_idx = set()
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if math.dist(sg[i][1], sg[j][1]) < ss_cut:
                ss_idx.add(sg[i][0])
                ss_idx.add(sg[j][0])
    breaks = set()
    for idx in range(1, cur):
        if idx in bbC and (idx + 1) in bbN:
            if math.dist(bbC[idx], bbN[idx + 1]) > 2.0:
                breaks.add(idx)
        else:
            breaks.add(idx)
    out, prevc = [], None
    for n, ln in enumerate(lines):
        if n in new_idx:
            idx = new_idx[n]
            if prevc is not None and idx != prevc and prevc in breaks:
                out.append("TER")
            prevc = idx
            rn = "CYX" if (ln[17:20] == "CYS" and idx in ss_idx) else ln[17:20]
            ln = ln[:17] + rn + ln[20:22] + f"{idx:4d}" + " " + ln[27:]
        out.append(ln)
    Path(out_pdb).write_text("\n".join(out) + "\nEND\n")
    return len(ss_idx) // 2


def fix_receptor(in_pdb: Path, out_pdb: Path) -> None:
    """Run PDBFixer to complete missing heavy atoms, replace nonstandard
    residues, and strip heterogens (waters/ions/ligands). Does not add
    hydrogens or model missing loops (see module docstring)."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(in_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # don't build in missing loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)


def prepare_receptor_pdb(raw_pdb: Path, out_pdb: Path, *, chain: str = "A",
                         lig_resname: str = "", drop_resseq_above: int | None = None,
                         tmp_dir: Path | None = None):
    """End-to-end: split one chain out of a raw (co-crystal) PDB, tidy it,
    and fix it with PDBFixer. Returns (out_pdb, ligand_lines, n_ss_bonds).

    `lig_resname` names the co-crystal ligand to extract alongside the
    receptor (used by `pocket.compute_box`); pass "" if there is none.
    """
    raw_pdb, out_pdb = Path(raw_pdb), Path(out_pdb)
    tmp_dir = Path(tmp_dir) if tmp_dir else out_pdb.parent
    tmp_dir.mkdir(parents=True, exist_ok=True)

    text = raw_pdb.read_text()
    receptor_lines, ligand_lines = split_receptor(text, chain, lig_resname, drop_resseq_above)
    if not receptor_lines:
        raise ValueError(f"{raw_pdb}: no ATOM lines found for chain {chain!r}")

    tidied = tmp_dir / f"{out_pdb.stem}_tidy.pdb"
    n_ss = tidy_receptor(receptor_lines, tidied)
    fix_receptor(tidied, out_pdb)
    return out_pdb, ligand_lines, n_ss
