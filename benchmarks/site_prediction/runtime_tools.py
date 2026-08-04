"""Runtime tool checks for site-prediction benchmarks."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKING_ROOT = REPO_ROOT / "docking"
DEFAULT_AUTODOCK_GPU_DIR = DOCKING_ROOT / "AUTODOCK_GPU_DIR"
DEFAULT_AUTOGRID4 = DEFAULT_AUTODOCK_GPU_DIR / "autogrid" / "autogrid4"
AUTODOCK_GPU_COMPILER = DOCKING_ROOT / "autodock-gpu-compiler.sh"


def first_command_token(command: str) -> str | None:
    parts = shlex.split(command)
    return parts[0] if parts else None


def is_executable_path(value: str | Path) -> bool:
    path = Path(value).expanduser()
    return path.is_file() and os.access(path, os.X_OK)


def command_is_available(command: str) -> bool:
    token = first_command_token(command)
    if not token:
        return False
    if "/" in token:
        return is_executable_path(token)
    return shutil.which(token) is not None


def _is_bundled_autogrid_path(path: Path) -> bool:
    try:
        return path.expanduser().resolve() == DEFAULT_AUTOGRID4.resolve()
    except FileNotFoundError:
        return path.expanduser().absolute() == DEFAULT_AUTOGRID4.absolute()


def build_bundled_autogrid4(*, numwi: int = 128) -> Path:
    """Compile/check the bundled AutoGrid source using Ultidock's normal build script."""
    if not AUTODOCK_GPU_COMPILER.is_file():
        raise FileNotFoundError(f"AutoDock-GPU compiler script not found: {AUTODOCK_GPU_COMPILER}")
    if not DEFAULT_AUTODOCK_GPU_DIR.is_dir():
        raise FileNotFoundError(f"AutoDock-GPU source directory not found: {DEFAULT_AUTODOCK_GPU_DIR}")

    env = os.environ.copy()
    env["DEVICE"] = "CPU"
    env["GPU_BACKEND"] = "CPU"
    env["GPU_TYPE"] = "CPU"
    env["NUMWI"] = str(numwi)
    print(f"[preflight] Building/checking bundled AutoGrid source: {DEFAULT_AUTOGRID4}")
    subprocess.run(
        ["bash", str(AUTODOCK_GPU_COMPILER), str(DEFAULT_AUTODOCK_GPU_DIR)],
        cwd=str(DEFAULT_AUTODOCK_GPU_DIR),
        env=env,
        check=True,
    )
    if not is_executable_path(DEFAULT_AUTOGRID4):
        raise RuntimeError(f"AutoGrid build finished but executable is missing: {DEFAULT_AUTOGRID4}")
    return DEFAULT_AUTOGRID4.resolve()


def resolve_autogrid4_bin(
    autogrid4_bin: str | Path,
    *,
    auto_build_bundled: bool = True,
    numwi: int = 128,
) -> str | None:
    """Return an executable AutoGrid path, building the bundled source if needed."""
    configured = Path(str(autogrid4_bin)).expanduser()
    if is_executable_path(configured):
        return str(configured.resolve())

    if _is_bundled_autogrid_path(configured):
        if auto_build_bundled:
            return str(build_bundled_autogrid4(numwi=numwi))
        return None

    token = str(autogrid4_bin)
    if "/" not in token:
        discovered = shutil.which(token)
        if discovered:
            return discovered
    return None
