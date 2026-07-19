"""Tests for GPU/CPU backend selection (dd_docking/gpu_backend.py).

These only exercise the pure-Python decision logic (`resolve_backend`,
`box_fits_gpu`) with `gpu_binary_available`/`platform.system` monkeypatched,
so they run everywhere regardless of whether Vina-GPU+ is actually built on
this machine. The real subprocess path (`dock_ligand_gpu`) is exercised
manually / in CI only where a built binary is present.
"""
import warnings
from pathlib import Path

import pytest

from dd_docking import gpu_backend


@pytest.fixture(autouse=True)
def _reset_warned(monkeypatch):
    monkeypatch.setattr(gpu_backend, "_warned_backends", set())


def test_box_fits_gpu():
    assert gpu_backend.box_fits_gpu([10.0, 20.0, 29.9])
    assert not gpu_backend.box_fits_gpu([10.0, 30.0, 10.0])
    assert not gpu_backend.box_fits_gpu([31.0, 5.0, 5.0])


def test_cpu_requested_always_cpu(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    assert gpu_backend.resolve_backend("cpu", [10, 10, 10]) == "cpu"


def test_auto_uses_gpu_when_available_and_box_fits(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    assert gpu_backend.resolve_backend("auto", [10, 10, 10]) == "gpu"


def test_auto_falls_back_to_cpu_without_binary(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: False)
    assert gpu_backend.resolve_backend("auto", [10, 10, 10]) == "cpu"


def test_auto_falls_back_to_cpu_when_box_too_large(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    assert gpu_backend.resolve_backend("auto", [10, 35, 10]) == "cpu"


def test_gpu_requested_without_binary_warns_and_falls_back(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: False)
    with pytest.warns(UserWarning, match="no Vina-GPU\\+ binary"):
        assert gpu_backend.resolve_backend("gpu", [10, 10, 10]) == "cpu"


def test_gpu_requested_with_oversized_box_warns_and_falls_back(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    with pytest.warns(UserWarning, match="cannot handle"):
        assert gpu_backend.resolve_backend("gpu", [10, 40, 10]) == "cpu"


def test_auto_falls_back_to_cpu_for_flexible_member(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    assert gpu_backend.resolve_backend("auto", [10, 10, 10], has_flex=True) == "cpu"


def test_gpu_requested_with_flex_warns_and_falls_back(monkeypatch):
    monkeypatch.setattr(gpu_backend, "gpu_binary_available", lambda: True)
    with pytest.warns(UserWarning, match="flexible side chains"):
        assert gpu_backend.resolve_backend("gpu", [10, 10, 10], has_flex=True) == "cpu"


def test_gpu_binary_unavailable_off_linux(monkeypatch):
    monkeypatch.setattr(gpu_backend.platform, "system", lambda: "Darwin")
    assert gpu_backend.gpu_binary_available() is False


def test_resolve_backend_rejects_unknown_value():
    with pytest.raises(ValueError):
        gpu_backend.resolve_backend("quantum", [10, 10, 10])


def test_dock_ligand_gpu_resolves_relative_receptor_paths(monkeypatch, tmp_path):
    # dock_ligand_gpu forces subprocess cwd to install_dir() (Kernel*_Opt.bin
    # are resolved relative to cwd, not the binary's own location), so a
    # relative rigid_pdbqt/flex_pdbqt -- as written into manifest.json by
    # every CLI example in the README -- must be made absolute first, or it
    # silently resolves against the wrong directory instead of the caller's.
    install = tmp_path / "install"
    install.mkdir()
    (install / gpu_backend.BIN_NAME).write_text("")
    monkeypatch.setattr(gpu_backend, "install_dir", lambda: install)
    monkeypatch.setattr(gpu_backend, "_raise_stack_limit", lambda: None)

    captured = {}

    def fake_run(args, **kwargs):
        captured["args"] = args
        raise FileNotFoundError  # short-circuit before actually exec'ing

    monkeypatch.setattr(gpu_backend.subprocess, "run", fake_run)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "rigid.pdbqt").write_text("")

    gpu_backend.dock_ligand_gpu("rigid.pdbqt", "LIGAND", [0, 0, 0], [10, 10, 10])

    receptor_arg = captured["args"][captured["args"].index("--receptor") + 1]
    assert Path(receptor_arg).is_absolute()
    assert receptor_arg == str(tmp_path / "rigid.pdbqt")


def test_warn_gpu_task_failed_warns_once():
    with pytest.warns(UserWarning, match="Vina-GPU\\+ docking failed"):
        gpu_backend.warn_gpu_task_failed("member1")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        gpu_backend.warn_gpu_task_failed("member1")  # second call: no warning
