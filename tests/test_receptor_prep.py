from dockinglib.receptor_prep import split_receptor, tidy_receptor


def _atom_line(serial, name, resname, chain, resnum, x, y, z, *, altloc=" ",
               icode=" ", element=None, rec="ATOM  "):
    """Build one fixed-column PDB ATOM/HETATM line (columns match the
    real PDB format used by `data/raw/*.pdb`)."""
    if element is None:
        element = name.strip()[0]
    name_field = f" {name:<3}" if len(name) < 4 else name
    return (f"{rec}{serial:>5} {name_field}{altloc}{resname:>3} {chain}{resnum:>4}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}")


def _ter_line(serial, resname, chain, resnum):
    return f"TER   {serial:>5}      {resname:>3} {chain}{resnum:>4}"


# --- split_receptor ---------------------------------------------------

def test_split_receptor_keeps_only_requested_chain():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "B", 1, 1.0, 0.0, 0.0),
    ])
    receptor, ligand = split_receptor(text, chain="A", lig_resname="LIG")
    assert len(receptor) == 1
    assert receptor[0][21] == "A"
    assert ligand == []


def test_split_receptor_extracts_matching_hetatm_ligand():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C1", "LIG", "A", 401, 5.0, 0.0, 0.0, rec="HETATM"),
        _atom_line(3, "O", "HOH", "A", 501, 9.0, 0.0, 0.0, rec="HETATM"),
    ])
    receptor, ligand = split_receptor(text, chain="A", lig_resname="LIG")
    assert len(receptor) == 1
    assert len(ligand) == 1
    assert ligand[0][17:20] == "LIG"


def test_split_receptor_keeps_primary_altloc_only():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, altloc="A"),
        _atom_line(2, "CA", "ALA", "A", 1, 0.1, 0.0, 0.0, altloc="B"),
        _atom_line(3, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, altloc=" "),
    ])
    receptor, _ = split_receptor(text, chain="A", lig_resname="")
    altlocs = [ln[16] for ln in receptor]
    assert altlocs == ["A", " "]


def test_split_receptor_drop_resseq_above():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "A", 50, 1.0, 0.0, 0.0),
        _atom_line(3, "CA", "VAL", "A", 200, 2.0, 0.0, 0.0),
    ])
    receptor, _ = split_receptor(text, chain="A", lig_resname="", drop_resseq_above=100)
    resnums = [int(ln[22:26]) for ln in receptor]
    assert resnums == [1, 50]


def test_split_receptor_keeps_ter_for_requested_chain():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _ter_line(2, "ALA", "A", 1),
        _atom_line(3, "CA", "GLY", "B", 1, 0.0, 0.0, 0.0),
        _ter_line(4, "GLY", "B", 1),
    ])
    receptor, _ = split_receptor(text, chain="A", lig_resname="")
    assert sum(1 for ln in receptor if ln.startswith("TER")) == 1


def test_split_receptor_ligand_restricted_to_first_matching_chain():
    text = "\n".join([
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C1", "LIG", "A", 401, 5.0, 0.0, 0.0, rec="HETATM"),
        _atom_line(3, "C1", "LIG", "B", 401, 5.0, 0.0, 0.0, rec="HETATM"),
    ])
    _, ligand = split_receptor(text, chain="A", lig_resname="LIG")
    assert len(ligand) == 1
    assert ligand[0][21] == "A"


# --- tidy_receptor ------------------------------------------------------

def test_tidy_receptor_renumbers_residues_sequentially(tmp_path):
    lines = [
        _atom_line(1, "CA", "ALA", "A", 5, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "A", 9, 1.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_receptor(lines, out_pdb)
    text = out_pdb.read_text()
    resnums = [int(ln[22:26]) for ln in text.splitlines() if ln[:6] in ("ATOM  ", "HETATM")]
    assert resnums == [1, 2]


def test_tidy_receptor_inserts_ter_at_backbone_break(tmp_path):
    # Residue 1 (N, C) and residue 2 (N, C) are far apart (> 2.0 A break).
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C", "ALA", "A", 1, 1.0, 0.0, 0.0),
        _atom_line(3, "N", "GLY", "A", 2, 50.0, 0.0, 0.0),
        _atom_line(4, "C", "GLY", "A", 2, 51.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_receptor(lines, out_pdb)
    out_lines = out_pdb.read_text().splitlines()
    assert "TER" in out_lines


def test_tidy_receptor_no_ter_when_backbone_is_continuous(tmp_path):
    # C(res1)-N(res2) close together (< 2.0 A): no chain break.
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C", "ALA", "A", 1, 1.3, 0.0, 0.0),
        _atom_line(3, "N", "GLY", "A", 2, 2.6, 0.0, 0.0),
        _atom_line(4, "C", "GLY", "A", 2, 3.9, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_receptor(lines, out_pdb)
    out_lines = out_pdb.read_text().splitlines()
    assert "TER" not in out_lines


def test_tidy_receptor_renames_disulfide_cys_to_cyx(tmp_path):
    # Two CYS residues with SG atoms close together (< default ss_cut=2.5 A)
    # should both be renamed to CYX.
    lines = [
        _atom_line(1, "SG", "CYS", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "SG", "CYS", "A", 2, 2.0, 0.0, 0.0),
        _atom_line(3, "SG", "CYS", "A", 3, 100.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    n_ss = tidy_receptor(lines, out_pdb)
    assert n_ss == 1
    out_lines = [ln for ln in out_pdb.read_text().splitlines() if ln[:6] in ("ATOM  ", "HETATM")]
    resnames = [ln[17:20] for ln in out_lines]
    # residues 1 and 2 (renumbered) are bonded -> CYX; residue 3 stays CYS.
    assert resnames == ["CYX", "CYX", "CYS"]


def test_tidy_receptor_returns_zero_disulfides_when_none_present(tmp_path):
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "A", 2, 1.3, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    n_ss = tidy_receptor(lines, out_pdb)
    assert n_ss == 0


def test_tidy_receptor_writes_end_record(tmp_path):
    lines = [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_receptor(lines, out_pdb)
    assert out_pdb.read_text().rstrip("\n").endswith("END")
