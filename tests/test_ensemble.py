import json

from dockinglib.ensemble import EnsembleMember, load_manifest, save_manifest


def _member(member_id="m1", flex_pdbqt="m1_flex.pdbqt"):
    return EnsembleMember(
        member_id=member_id,
        receptor_pdb=f"{member_id}_fixed.pdb",
        rigid_pdbqt=f"{member_id}_rigid.pdbqt",
        flex_pdbqt=flex_pdbqt,
        center=[1.0, 2.0, 3.0],
        size=[20.0, 21.0, 22.0],
        flexres=["A:10", "A:20"],
    )


def test_save_manifest_writes_json_list(tmp_path):
    out = tmp_path / "manifest.json"
    save_manifest(out, [_member("m1"), _member("m2")])

    data = json.loads(out.read_text())
    assert isinstance(data, list)
    assert len(data) == 2
    assert data[0]["member_id"] == "m1"
    assert data[0]["flexres"] == ["A:10", "A:20"]


def test_load_manifest_roundtrip(tmp_path):
    out = tmp_path / "manifest.json"
    original = [_member("m1"), _member("m2", flex_pdbqt=None)]
    save_manifest(out, original)

    loaded = load_manifest(out)

    assert loaded == original
    assert all(isinstance(m, EnsembleMember) for m in loaded)


def test_load_manifest_preserves_none_flex_pdbqt(tmp_path):
    out = tmp_path / "manifest.json"
    save_manifest(out, [_member("rigid_only", flex_pdbqt=None)])

    loaded = load_manifest(out)

    assert loaded[0].flex_pdbqt is None


def test_load_manifest_reads_real_repo_fixture():
    # Reuse the manifest already checked into data/ensemble/ from an
    # earlier real dockinglib-prep run, to make sure load_manifest can
    # parse a real (not just synthetic) manifest.json.
    from pathlib import Path

    manifest_path = Path(__file__).resolve().parent.parent / "data" / "ensemble" / "manifest.json"
    if not manifest_path.exists():
        import pytest
        pytest.skip("data/ensemble/manifest.json not present in this checkout")

    members = load_manifest(manifest_path)
    assert len(members) >= 1
    assert all(isinstance(m, EnsembleMember) for m in members)
    assert all(m.member_id for m in members)
