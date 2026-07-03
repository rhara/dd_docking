from dockinglib.screening import EnsembleHit, rank_key, write_results_csv


def _hit(ligand_id, best_affinity, best_member="m1"):
    return EnsembleHit(
        ligand_id=ligand_id, smiles="C", best_member=best_member,
        best_affinity=best_affinity, member_affinities={best_member: best_affinity},
    )


def test_rank_key_sorts_lowest_affinity_first():
    hits = [_hit("b", -5.0), _hit("a", -9.0), _hit("c", -1.0)]
    hits.sort(key=rank_key)
    assert [h.ligand_id for h in hits] == ["a", "b", "c"]


def test_rank_key_sorts_failed_ligands_last():
    hits = [_hit("ok", -7.0), _hit("failed", None)]
    hits.sort(key=rank_key)
    assert [h.ligand_id for h in hits] == ["ok", "failed"]


def test_write_results_csv_columns(tmp_path):
    import pandas as pd

    hits = [_hit("lig1", -8.5, "m1"), _hit("lig2", -6.0, "m2")]
    out_csv = tmp_path / "ranked.csv"
    write_results_csv(str(out_csv), hits)

    df = pd.read_csv(out_csv)
    assert list(df["ligand_id"]) == ["lig1", "lig2"]
    assert list(df["rank"]) == [1, 2]
    assert list(df["best_affinity"]) == [-8.5, -6.0]
    assert "affinity[m1]" in df.columns
