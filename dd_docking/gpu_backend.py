"""Optional GPU docking backend: Vina-GPU+ (DeltaGroupNJUPT/Vina-GPU-2.0),
built and installed by `scripts/build_vina_gpu.sh`.

Linux-only, and rigid-receptor-only: Vina-GPU+ has no support at all for
flexible side chains (see `resolve_backend`'s docstring). `resolve_backend`
is the single place that decides whether a given docking task can/should
run on the GPU; everything else (macOS, Windows, no binary installed, a
member with flexible residues, or a docking box too large for the OpenCL
kernel's fixed limit) transparently falls back to the CPU QuickVina2 backend
in `docking.py`, so callers never need their own platform checks.
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Optional, Sequence, Tuple

from . import io

BIN_NAME = "Vina-GPU"
KERNEL_FILES = ("Kernel1_Opt.bin", "Kernel2_Opt.bin")

# Fixed limitation of the Vina-GPU+ OpenCL kernel (see upstream README):
# search box must be under 30 Angstrom in each dimension.
MAX_BOX_DIM = 30.0

# 8 MiB, per upstream README ("At least 8M stack size is needed" on Linux).
_REQUIRED_STACK_BYTES = 8 * 1024 * 1024

_warned_backends: set = set()


def install_dir() -> Path:
    """Directory expected to hold the `Vina-GPU` binary and its two
    Kernel*_Opt.bin files, installed by `scripts/build_vina_gpu.sh`.
    """
    env_dir = os.environ.get("DD_DOCKING_VINA_GPU_DIR")
    if env_dir:
        return Path(env_dir)
    return Path(sys.prefix) / "share" / "dd_docking" / "vina-gpu"


def gpu_binary_available() -> bool:
    """Whether a Vina-GPU+ binary (with its kernel binaries) is installed
    and usable on this platform. Always False off Linux.
    """
    if platform.system() != "Linux":
        return False
    d = install_dir()
    return all((d / name).exists() for name in (BIN_NAME, *KERNEL_FILES))


def box_fits_gpu(size: Sequence[float]) -> bool:
    return all(dim < MAX_BOX_DIM for dim in size)


def _warn_once(key: str, message: str) -> None:
    if key not in _warned_backends:
        _warned_backends.add(key)
        warnings.warn(message, stacklevel=3)


def resolve_backend(requested: str, size: Sequence[float], *, has_flex: bool = False) -> str:
    """Decide "cpu" or "gpu" for a docking box of the given `size`.

    `requested` is "cpu" (always CPU), "gpu" (GPU if at all possible, else
    warn once and fall back to CPU), or "auto" (GPU if available and the
    box fits the kernel's size limit, else CPU silently).

    `has_flex` must be True if this member has any flexible side chains --
    Vina-GPU+'s `main_procedure_cl.cpp` asserts `m.num_other_pairs() == 0`,
    which is nonzero for any ligand-flex/flex-flex/flex-inflex interaction,
    so it only ever supports fully rigid receptors. This is a hard limit of
    the tool itself (not a box-size or driver issue), so it's checked first
    and unconditionally routes to CPU.
    """
    if requested not in ("auto", "cpu", "gpu"):
        raise ValueError(f"resolve_backend: unknown backend {requested!r}")

    if requested == "cpu":
        return "cpu"

    if has_flex:
        if requested == "gpu":
            _warn_once(
                "gpu-no-flex-support",
                "dd_docking: --backend gpu requested but this ensemble member has "
                "flexible side chains -- Vina-GPU+ only supports rigid receptors "
                "(it asserts num_other_pairs() == 0) -- falling back to CPU (QuickVina2).",
            )
        return "cpu"

    if not gpu_binary_available():
        if requested == "gpu":
            _warn_once(
                "gpu-unavailable",
                "dd_docking: --backend gpu requested but no Vina-GPU+ binary found "
                f"for this platform/{install_dir()} -- falling back to CPU (QuickVina2). "
                "Run scripts/build_vina_gpu.sh on Linux to build it.",
            )
        return "cpu"

    if not box_fits_gpu(size):
        if requested == "gpu":
            _warn_once(
                "gpu-box-too-large",
                f"dd_docking: docking box {list(size)} has a dimension >= {MAX_BOX_DIM} A, "
                "which Vina-GPU+'s OpenCL kernel cannot handle -- falling back to CPU "
                "(QuickVina2) for this ensemble member.",
            )
        return "cpu"

    return "gpu"


def warn_gpu_task_failed(member_id: str) -> None:
    """Called by screening.py when `dock_ligand_gpu` returns None for a
    member `resolve_backend` said should use the GPU -- warns once per
    process, since this usually means the binary/driver/OpenCL kernel is
    broken on this machine (some GPU generations hit upstream Vina-GPU+
    kernel build failures) rather than a one-off bad ligand.
    """
    _warn_once(
        "gpu-task-failed",
        f"dd_docking: Vina-GPU+ docking failed for ensemble member {member_id!r} "
        "(binary/driver/OpenCL kernel issue) -- falling back to CPU (QuickVina2) for this "
        "and any further failures this run.",
    )


def dock_ligand_gpu(
    rigid_pdbqt: str,
    ligand_pdbqt: str,
    center: Sequence[float],
    size: Sequence[float],
    *,
    flex_pdbqt: Optional[str] = None,
    seed: int = 0,
    n_poses: int = 5,
    thread: int = 8000,
    search_depth: Optional[int] = None,
) -> Optional[Tuple[float, str]]:
    """Dock one ligand with Vina-GPU+. Mirrors `docking.dock_ligand`'s
    return contract: (best_affinity_kcal_mol, poses_pdbqt) or None on
    failure.
    """
    d = install_dir()
    binary = d / BIN_NAME
    # The subprocess's cwd is forced to `d` below (Kernel*_Opt.bin are
    # resolved relative to cwd, not the binary's location), so a relative
    # rigid_pdbqt/flex_pdbqt from the caller -- e.g. manifest.json entries
    # are written relative to the ensemble output dir in every CLI example
    # in the README -- would otherwise resolve against the wrong directory
    # and silently fail to open.
    rigid_pdbqt = str(Path(rigid_pdbqt).resolve())
    if flex_pdbqt:
        flex_pdbqt = str(Path(flex_pdbqt).resolve())

    with tempfile.TemporaryDirectory(prefix="dd_docking_gpu_") as tmp:
        tmp_path = Path(tmp)
        ligand_path = tmp_path / "ligand.pdbqt"
        ligand_path.write_text(ligand_pdbqt)
        out_path = tmp_path / "out.pdbqt"

        args = [
            str(binary),
            "--receptor", str(rigid_pdbqt),
            "--ligand", str(ligand_path),
            "--out", str(out_path),
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--thread", str(thread),
            "--num_modes", str(n_poses),
            "--seed", str(seed),
        ]
        if flex_pdbqt:
            args += ["--flex", str(flex_pdbqt)]
        if search_depth is not None:
            args += ["--search_depth", str(search_depth)]

        env = dict(os.environ)
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            filter(None, [str(d), env.get("LD_LIBRARY_PATH")])
        )
        try:
            subprocess.run(
                args, cwd=d, env=env, capture_output=True, text=True, check=True,
                preexec_fn=_raise_stack_limit,
            )
        except (subprocess.CalledProcessError, OSError):
            return None

        if not out_path.exists():
            return None
        poses_pdbqt = out_path.read_text()

    affinity = io.parse_vina_affinity(poses_pdbqt)
    if affinity is None:
        return None
    return affinity, poses_pdbqt


def _raise_stack_limit() -> None:
    """Vina-GPU+ needs an 8 MiB stack on Linux (upstream README); raise the
    child process's limit before exec so callers don't need `ulimit -s`.
    """
    import resource

    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    target = _REQUIRED_STACK_BYTES
    if hard != resource.RLIM_INFINITY:
        target = min(target, hard)
    if soft < target:
        resource.setrlimit(resource.RLIMIT_STACK, (target, hard))
