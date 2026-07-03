from pathlib import Path

from dockinglib.pocket import Residue, compute_box, find_pocket_residues, format_flexres

# Two ATOM lines (chain A, resnum 10 and 20) and one far-away line (resnum 99).
_PDB_TEXT = (
    "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      2  CA  GLY A  20       2.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      3  CA  VAL A  99      50.000   0.000   0.000  1.00  0.00           C\n"
)

_LIGAND_LINES = [
    "HETATM    1  C1  LIG A 401       0.500   0.000   0.000  1.00  0.00           C",
    "HETATM    2  C2  LIG A 401       1.500   1.000  -1.000  1.00  0.00           C",
]


def test_compute_box_centers_on_ligand_with_padding():
    center, size = compute_box(_LIGAND_LINES, padding=2.0)
    assert center == [1.0, 0.5, -0.5]
    assert size == [1.0 + 4.0, 1.0 + 4.0, 1.0 + 4.0]


def test_find_pocket_residues_within_cutoff(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(_PDB_TEXT)

    ref_coords = [(0.0, 0.0, 0.0)]
    residues = find_pocket_residues(pdb_path, ref_coords, cutoff=3.0)

    assert residues == [Residue("A", 10), Residue("A", 20)]


def test_find_pocket_residues_caps_to_closest(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(_PDB_TEXT)

    ref_coords = [(0.0, 0.0, 0.0)]
    # Both residues are within cutoff, but max_residues=1 should keep only
    # the closer one (resnum 10, distance 0.0) over resnum 20 (distance 2.0).
    residues = find_pocket_residues(pdb_path, ref_coords, cutoff=3.0, max_residues=1)

    assert residues == [Residue("A", 10)]


def test_format_flexres_is_explicit_chain_prefixed():
    residues = [Residue("A", 10), Residue("A", 20), Residue("B", 5)]
    assert format_flexres(residues) == "A:10,A:20,B:5"


def test_format_flexres_empty():
    assert format_flexres([]) == ""
