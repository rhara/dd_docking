from dd_docking.io import read_results_csv, write_results_csv
from dd_docking.ligand_prep import read_smi


def test_write_read_results_csv_roundtrip(tmp_path):
    rows = [
        {"rank": 1, "ligand_id": "lig1", "best_affinity": -8.5},
        {"rank": 2, "ligand_id": "lig2", "best_affinity": -6.0},
    ]
    out_csv = tmp_path / "results.csv"
    write_results_csv(str(out_csv), rows)

    df = read_results_csv(str(out_csv))
    assert list(df["ligand_id"]) == ["lig1", "lig2"]
    assert list(df["best_affinity"]) == [-8.5, -6.0]


def test_read_smi_parses_id_and_skips_comments(tmp_path):
    smi_path = tmp_path / "ligands.smi"
    smi_path.write_text("# comment\nCCO\tethanol\n\nc1ccccc1  benzene\n")

    ligands = read_smi(smi_path)

    assert [(l.ligand_id, l.smiles) for l in ligands] == [
        ("ethanol", "CCO"),
        ("benzene", "c1ccccc1"),
    ]


def test_read_smi_rejects_missing_id(tmp_path):
    import pytest

    smi_path = tmp_path / "bad.smi"
    smi_path.write_text("CCO\n")

    with pytest.raises(ValueError):
        read_smi(smi_path)
