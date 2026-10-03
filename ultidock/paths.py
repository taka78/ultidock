"""Keep installed packages read-only and managed workflows in a writable home."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from zipfile import ZipFile

from ultidock import __version__


def source_root() -> Path | None:
    root = Path(__file__).resolve().parents[1]
    if (root / "pyproject.toml").is_file() and (root / "docking/run.py").is_file():
        return root
    return None


def workspace_root() -> Path:
    """Use an explicit home, the development checkout, or the user data directory."""
    if os.environ.get("ULTIDOCK_HOME"):
        return Path(os.environ["ULTIDOCK_HOME"]).expanduser().resolve()
    checkout = source_root()
    if checkout is not None:
        return checkout
    data_home = Path(os.environ.get("XDG_DATA_HOME", "~/.local/share")).expanduser()
    return (data_home / "ultidock").resolve()


def _link_native_programs(root: Path, checkout: Path | None) -> None:
    """Adapt PATH programs to the directory layout expected by legacy scripts."""
    programs = {
        "vina": "docking/VINA_DIR/bin/vina",
        "vina_split": "docking/VINA_DIR/bin/vina_split",
        "autogrid4": "docking/AUTODOCK_GPU_DIR/autogrid/autogrid4",
    }
    for name, relative in programs.items():
        target = root / relative
        if target.exists() or target.is_symlink():
            continue
        bundled = checkout / relative if checkout else None
        program = str(bundled) if bundled and os.access(bundled, os.X_OK) else shutil.which(name)
        if program:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(Path(program).resolve())


def application_root() -> Path:
    """Materialize bundled scripts/assets once; never copy configs or run outputs."""
    root = workspace_root()
    checkout = source_root()
    if root == checkout:
        return root
    package_root = Path(__file__).resolve().parents[1]
    # A subprocess executing a materialized workflow already has its resources.
    if root == package_root and (root / ".ultidock-runtime.json").is_file():
        return root
    archive_path = Path(__file__).with_name("runtime.zip")
    if checkout is None and not archive_path.is_file():
        raise RuntimeError(
            "Ultidock runtime resources are missing; reinstall the ultidock package."
        )
    root.mkdir(parents=True, exist_ok=True)
    # Workflows and their native build scripts currently require a POSIX host.
    import fcntl

    with (root / ".ultidock-runtime.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if checkout is not None:
            from ultidock.runtime_files import runtime_files

            resources = [
                (path.relative_to(checkout), path.read_bytes(), path.stat().st_mode & 0o777)
                for path in runtime_files(checkout)
            ]
        else:
            with ZipFile(archive_path) as archive:
                resources = [
                    (Path(info.filename), archive.read(info), (info.external_attr >> 16) & 0o777)
                    for info in archive.infolist()
                    if not info.is_dir()
                ]
        digest = hashlib.sha256()
        for relative, content, mode in resources:
            digest.update(relative.as_posix().encode())
            digest.update(content)
        revision = digest.hexdigest()
        marker = root / ".ultidock-runtime.json"
        try:
            current = json.loads(marker.read_text())
        except (FileNotFoundError, ValueError):
            current = {}
        if current.get("revision") != revision:
            for relative, content, mode in resources:
                if relative.is_absolute() or ".." in relative.parts:
                    raise RuntimeError(f"Invalid bundled resource: {relative}")
                target = root / relative
                if root not in target.resolve().parents:
                    raise RuntimeError(f"Resource path escapes ULTIDOCK_HOME: {target}")
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.is_file() and target.read_bytes() == content:
                    continue
                with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as temporary:
                    temporary.write(content)
                    temporary_path = Path(temporary.name)
                temporary_path.chmod(mode or 0o644)
                temporary_path.replace(target)
            marker.write_text(json.dumps({"version": __version__, "revision": revision}) + "\n")
        _link_native_programs(root, checkout)
    return root
