"""Exercise a built release outside the checkout, with an independent workspace.

The temporary venv shares the test environment's Python dependencies; building
and installing the package itself never requires a network connection.
"""

from __future__ import annotations

import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

from ultidock.runtime_files import runtime_files


def _run(command, *, cwd, env=None):
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def test_sdist_wheel_installs_and_runs_outside_checkout(tmp_path):
    repo = Path(__file__).resolve().parents[2]
    source = tmp_path / "source"
    source.mkdir()
    for path in [repo / "pyproject.toml", *runtime_files(repo)]:
        target = source / path.relative_to(repo)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    # These local files must never enter either a release or a fresh workspace.
    (source / "docking/config.py").write_text("DO_NOT_SHIP = True\n")
    results = source / "examples/sert-escitalopram/workspace/old/results.csv"
    results.parent.mkdir(parents=True)
    results.write_text("DO_NOT_SHIP\n")

    dist = tmp_path / "dist"
    dist.mkdir()
    build_env = dict(os.environ)
    build_env.pop("PYTHONPATH", None)
    _run(
        [
            sys.executable,
            "-c",
            f"from setuptools.build_meta import build_sdist; build_sdist({str(dist)!r})",
        ],
        cwd=source,
        env=build_env,
    )
    sdist = next(dist.glob("ultidock-*.tar.gz"))
    with tarfile.open(sdist) as archive:
        names = archive.getnames()
        assert any(name.endswith("ad4_shared/paramdat2h.csh") for name in names)
        assert not any(
            name.endswith("docking/config.py") or "/workspace/" in name for name in names
        )

    _run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(sdist),
            "--no-deps",
            "--no-build-isolation",
            "--wheel-dir",
            str(dist),
        ],
        cwd=tmp_path,
        env=build_env,
    )
    wheel = next(dist.glob("ultidock-*.whl"))
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        assert "docking/pocket_boxes.py" in names
        assert "benchmarks/baseline/common.py" in names
        assert "docking/config.py" not in names
        assert not any(name.startswith("molguard/tests/") for name in names)
        with zipfile.ZipFile(io.BytesIO(archive.read("ultidock/runtime.zip"))) as runtime:
            assert "examples/quickstart/data/receptor.pdb" in runtime.namelist()
            assert (
                "docking/AUTODOCK_GPU_DIR/autogrid/ad4_shared/paramdat2h.csh" in runtime.namelist()
            )
            assert not any(
                "/workspace/" in name or name.endswith("config.py") for name in runtime.namelist()
            )

    venv = tmp_path / "venv"
    _run([sys.executable, "-m", "venv", "--system-site-packages", str(venv)], cwd=tmp_path)
    python = str(venv / "bin/python")
    _run([python, "-m", "pip", "install", "--no-deps", str(wheel)], cwd=tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    home = tmp_path / "managed_home"
    env = {**build_env, "ULTIDOCK_HOME": str(home), "PYTHONDONTWRITEBYTECODE": "1"}
    # Ensure imports really come from this wheel, rather than the test checkout.
    origins = json.loads(
        _run(
            [
                python,
                "-c",
                "import json, ultidock, docking.pocket_boxes, benchmarks.baseline.common; "
                "from importlib.metadata import version; "
                "print(json.dumps([version('ultidock'), ultidock.__file__, "
                "docking.pocket_boxes.__file__, benchmarks.baseline.common.__file__]))",
            ],
            cwd=outside,
            env=env,
        )
    )
    assert origins[0] == "1.1.1"
    assert all(str(venv) in origin for origin in origins[1:])
    installed = Path(origins[1]).parents[1]
    initial_files = {path.relative_to(installed) for path in installed.rglob("*") if path.is_file()}
    ultidock = str(venv / "bin/ultidock")
    assert "Commands:" in _run([ultidock, "--help"], cwd=outside, env=env)
    assert not home.exists()
    assert "version 1.1.1" in _run([python, "-m", "ultidock", "--version"], cwd=outside, env=env)
    assert "version 1.1.1" in _run([str(venv / "bin/molguard"), "--version"], cwd=outside, env=env)
    assert "sert-escitalopram" in _run([ultidock, "example", "list"], cwd=outside, env=env)
    _run([ultidock, "example", "run", "quickstart"], cwd=outside, env=env)
    quickstart = home / "examples/quickstart/workspace/quickstart_run"
    assert (quickstart / "report.html").is_file()
    _run([ultidock, "report", str(quickstart)], cwd=outside, env=env)
    for mode in ("known-site", "cavity", "blind"):
        output = home / "runs" / mode
        args = [ultidock, mode, "--output-dir", str(output), "--dry-run"]
        if mode == "known-site":
            args += ["--center", "1,2,3"]
        _run(args, cwd=outside, env=env)
        assert (output / "run_config.yaml").is_file()
    _run([ultidock, "clean"], cwd=outside, env=env)
    assert (quickstart / "report.html").is_file()
    assert initial_files == {
        path.relative_to(installed) for path in installed.rglob("*") if path.is_file()
    }

    # Editable development installs must also expose both command entry points.
    _run(
        [python, "-m", "pip", "install", "--no-deps", "--no-build-isolation", "-e", str(source)],
        cwd=outside,
    )
    _run([ultidock, "--help"], cwd=outside, env=env)
    _run([str(venv / "bin/molguard"), "--version"], cwd=outside, env=env)
