import numpy as np
import pytest

from dockinglib.refine_md import _summarize_rmsd


def test_summarize_rmsd_stable_when_low_and_flat():
    rmsd = np.array([0.5, 0.8, 1.0, 0.9, 1.1])
    summary = _summarize_rmsd(rmsd)

    assert summary["n_frames"] == 5
    assert summary["rmsd_mean"] == pytest.approx(0.86, abs=0.01)
    assert summary["rmsd_final"] == pytest.approx(1.1, abs=1e-9)
    assert summary["rmsd_max"] == pytest.approx(1.1, abs=1e-9)
    assert summary["stable"] is True


def test_summarize_rmsd_unstable_when_mean_too_high():
    # Mean RMSD (5.0) exceeds the default 3.0 A cutoff.
    rmsd = np.array([5.0, 5.0, 5.0])
    summary = _summarize_rmsd(rmsd)
    assert summary["stable"] is False


def test_summarize_rmsd_unstable_when_final_frame_drifted():
    # Mean is low but the final frame drifted past the 4.0 A cutoff.
    rmsd = np.array([0.5, 0.5, 0.5, 0.5, 6.0])
    summary = _summarize_rmsd(rmsd)
    assert summary["rmsd_mean"] < 3.0
    assert summary["stable"] is False


def test_summarize_rmsd_boundary_is_strict_less_than():
    # Exactly at the cutoffs should NOT count as stable (strict <).
    rmsd = np.array([3.0, 3.0, 4.0])
    summary = _summarize_rmsd(rmsd)
    assert summary["rmsd_mean"] == pytest.approx(10.0 / 3.0, abs=0.01)
    assert summary["stable"] is False


def test_summarize_rmsd_respects_custom_cutoffs():
    rmsd = np.array([5.0, 5.0, 5.0])
    summary = _summarize_rmsd(rmsd, stable_mean_cutoff=10.0, stable_final_cutoff=10.0)
    assert summary["stable"] is True


def test_summarize_rmsd_single_frame():
    rmsd = np.array([1.5])
    summary = _summarize_rmsd(rmsd)
    assert summary["n_frames"] == 1
    assert summary["rmsd_mean"] == summary["rmsd_final"] == summary["rmsd_max"] == 1.5


# --- heavy, environment-dependent paths: skipped ------------------------

@pytest.mark.skip(
    reason="requires a real multi-minute OpenMM MD run (minimize + equilibrate "
           "+ production) with a real protonated receptor PDB and docked pose "
           "PDBQT; not practical to run in the unit test suite."
)
def test_refine_one_hit_full_md_run():
    pass


@pytest.mark.skip(
    reason="end-to-end refine_top_hits() drives refine_one_hit() per row and "
           "therefore needs the same real multi-minute MD run plus a real "
           "ranked_results.csv referencing existing receptor/pose files; not "
           "practical to run in the unit test suite."
)
def test_refine_top_hits_full_pipeline():
    pass
