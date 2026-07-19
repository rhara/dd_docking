"""Tests for the QuickVina2 wrapper in docking.py.

These reuse the real rigid/flex PDBQT receptor fixtures already checked
into data/ensemble/ (produced by an earlier real dd_docking-prep run) so
we exercise the actual `qvina2` binary end-to-end, but keep exhaustiveness
low and the ligand tiny so the whole suite still runs in well under a
minute.
"""
import shutil
from pathlib import Path

import pytest

from dd_docking.docking import dock_ligand

_DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "ensemble"
_RIGID_PDBQT = _DATA_DIR / "3ert_rigid.pdbqt"
_FLEX_PDBQT = _DATA_DIR / "3ert_flex.pdbqt"
_CENTER = [30.282, -1.913, 24.206]
_SIZE = [28.448, 24.756, 28.687]

pytestmark = [
    pytest.mark.skipif(shutil.which("qvina2") is None, reason="qvina2 binary not on PATH"),
    pytest.mark.skipif(
        not _RIGID_PDBQT.exists(),
        reason="data/ensemble/3ert_rigid.pdbqt fixture not present in this checkout",
    ),
]


def test_dock_ligand_returns_none_on_invalid_pdbqt():
    result = dock_ligand(
        str(_RIGID_PDBQT), "not a valid pdbqt string", _CENTER, _SIZE,
        seed=0, exhaustiveness=1, n_poses=1, cpu=1,
    )
    assert result is None


def test_dock_ligand_real_small_ligand_returns_affinity_and_poses():
    pytest.importorskip("meeko")
    from dd_docking.ligand_prep import prepare_ligand_pdbqt

    pdbqt = prepare_ligand_pdbqt("CCO", seed=0)  # ethanol: tiny, embeds fast
    assert pdbqt is not None

    result = dock_ligand(
        str(_RIGID_PDBQT), pdbqt, _CENTER, _SIZE,
        seed=0, exhaustiveness=1, n_poses=1, cpu=1,
    )

    assert result is not None
    affinity, poses_pdbqt = result
    assert isinstance(affinity, float)
    assert isinstance(poses_pdbqt, str)
    assert "ATOM" in poses_pdbqt or "HETATM" in poses_pdbqt


def test_dock_ligand_with_flex_pdbqt_returns_affinity():
    pytest.importorskip("meeko")
    from dd_docking.ligand_prep import prepare_ligand_pdbqt

    pdbqt = prepare_ligand_pdbqt("CCO", seed=0)
    assert pdbqt is not None

    result = dock_ligand(
        str(_RIGID_PDBQT), pdbqt, _CENTER, _SIZE,
        flex_pdbqt=str(_FLEX_PDBQT), seed=0, exhaustiveness=1, n_poses=1, cpu=1,
    )

    assert result is not None
    affinity, poses_pdbqt = result
    assert isinstance(affinity, float)
