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


def regularize_carboxylate_geometry(topology, positions, *, bad_distance_nm: float = 0.20):
    """Fix a PDBFixer `addMissingAtoms()` quirk: a freshly-added carboxylate
    partner oxygen (backbone OXT at a chain terminus, or Asp/Glu OD2/OE2)
    can end up at a chemically impossible angle from its sibling oxygen --
    observed as close as ~0.10-0.19 nm apart instead of the ~0.22 nm a
    symmetric carboxylate requires, with individually normal C-O bond
    lengths (the angle is wrong, not the bond length). A real, well-formed
    carboxylate's two oxygens are consistently ~0.21-0.23 nm apart, so
    0.20 nm cleanly separates genuine defects from normal geometry. Meeko's
    per-residue, purely-distance-based bond perception (`mk_prepare_receptor.py`,
    via RDKit's `DetermineConnectivity`) then misreads the short O-O
    distance as a real bond and crashes with an oxygen valence error.

    Detects any carbon bonded to exactly two oxygens closer to each other
    than `bad_distance_nm`. The carbon is sp2 (three substituents: the two
    oxygens plus one "stem" neighbor -- CA for a backbone terminus, CG/CB
    for a Glu/Asp side chain), so the O-C-O bisector is fixed by trigonal
    geometry: opposite the stem bond, independent of the (broken) original
    oxygen positions. Rebuilds both oxygens there, ~124 degrees apart,
    preserving each one's own original C-O bond length. The rotation
    around the stem axis is not chemically constrained by this alone and
    is picked arbitrarily but consistently; that's fine here since the
    goal is just a clash-free, correctly-angled geometry for Meeko/Vina,
    not reproducing an unobserved crystallographic OXT position exactly.
    Returns a new positions array; does not mutate `positions` in place.
    """
    import numpy as np
    from openmm import unit

    pos = np.array(positions.value_in_unit(unit.nanometer))
    atoms = list(topology.atoms())
    neighbors = {a.index: set() for a in atoms}
    for a1, a2 in topology.bonds():
        neighbors[a1.index].add(a2.index)
        neighbors[a2.index].add(a1.index)

    half_angle = math.radians(62.0)
    n_fixed = 0
    for atom in atoms:
        if atom.element is None or atom.element.symbol != "C":
            continue
        o_idx = [i for i in neighbors[atom.index] if atoms[i].element and atoms[i].element.symbol == "O"]
        if len(o_idx) != 2:
            continue
        i1, i2 = o_idx
        c = pos[atom.index]
        if np.linalg.norm(pos[i1] - pos[i2]) >= bad_distance_nm:
            continue
        len1, len2 = np.linalg.norm(pos[i1] - c), np.linalg.norm(pos[i2] - c)

        other = [i for i in neighbors[atom.index] if i not in (i1, i2)]
        if not other:
            continue  # no stem bond to derive the plane from; leave as-is
        stem = pos[other[0]] - c
        stem_len = np.linalg.norm(stem)
        if stem_len < 1e-6:
            continue
        axis = -stem / stem_len  # sp2 trigonal geometry: O-C-O bisector is opposite the stem bond

        helper = pos[other[1]] - c if len(other) > 1 else np.array([1.0, 0.0, 0.0])
        perp = helper - np.dot(helper, axis) * axis
        if np.linalg.norm(perp) < 1e-6:
            helper = np.array([0.0, 1.0, 0.0])
            perp = helper - np.dot(helper, axis) * axis
        perp = perp / np.linalg.norm(perp)

        pos[i1] = c + len1 * (math.cos(half_angle) * axis + math.sin(half_angle) * perp)
        pos[i2] = c + len2 * (math.cos(half_angle) * axis - math.sin(half_angle) * perp)
        n_fixed += 1

    if n_fixed:
        print(f"[dd_docking] regularized {n_fixed} carboxylate geometry defect(s) from PDBFixer", flush=True)
    return unit.Quantity(pos, unit.nanometer)


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
    positions = regularize_carboxylate_geometry(fixer.topology, fixer.positions)
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, positions, fh, keepIds=True)


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
