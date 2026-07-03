"""Tests for the Vina wrapper in docking.py.

These reuse the real rigid/flex PDBQT receptor fixtures already checked
into data/ensemble/ (produced by an earlier real dd_docking-prep run) so
we exercise the actual `vina` package end-to-end, but keep exhaustiveness
low and the ligand tiny so the whole suite still runs in well under a
minute.
"""
from pathlib import Path

import pytest

vina = pytest.importorskip("vina")

from dd_docking.docking import dock_ligand, make_vina  # noqa: E402

_DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "ensemble"
_RIGID_PDBQT = _DATA_DIR / "6w63_rigid.pdbqt"
_FLEX_PDBQT = _DATA_DIR / "6w63_flex.pdbqt"
_CENTER = [-19.343, 18.299, -27.242]
_SIZE = [17.944, 21.143, 22.978]

pytestmark = pytest.mark.skipif(
    not _RIGID_PDBQT.exists(),
    reason="data/ensemble/6w63_rigid.pdbqt fixture not present in this checkout",
)


def test_make_vina_builds_rigid_only_instance():
    v = make_vina(str(_RIGID_PDBQT), _CENTER, _SIZE, seed=0, exhaustiveness=4, cpu=1)
    assert v._exhaustiveness == 4


def test_make_vina_builds_with_flex_pdbqt():
    v = make_vina(
        str(_RIGID_PDBQT), _CENTER, _SIZE,
        flex_pdbqt=str(_FLEX_PDBQT), seed=0, exhaustiveness=4, cpu=1,
    )
    assert v._exhaustiveness == 4


def test_dock_ligand_returns_none_on_invalid_pdbqt():
    v = make_vina(str(_RIGID_PDBQT), _CENTER, _SIZE, seed=0, exhaustiveness=1, cpu=1)
    result = dock_ligand(v, "not a valid pdbqt string", n_poses=1)
    assert result is None


def test_dock_ligand_real_small_ligand_returns_affinity_and_poses():
    ligand_prep = pytest.importorskip("meeko")
    from dd_docking.ligand_prep import prepare_ligand_pdbqt

    pdbqt = prepare_ligand_pdbqt("CCO", seed=0)  # ethanol: tiny, embeds fast
    assert pdbqt is not None

    v = make_vina(str(_RIGID_PDBQT), _CENTER, _SIZE, seed=0, exhaustiveness=1, cpu=1)
    result = dock_ligand(v, pdbqt, n_poses=1)

    assert result is not None
    affinity, poses_pdbqt = result
    assert isinstance(affinity, float)
    assert isinstance(poses_pdbqt, str)
    assert "ATOM" in poses_pdbqt or "HETATM" in poses_pdbqt
