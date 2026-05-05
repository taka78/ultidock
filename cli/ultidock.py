"""Ultidock command-line interface for workflow orchestration."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import click

from cli.readme import readme_hint
from molguard import __version__

FORWARD_CONTEXT = {"ignore_unknown_options": True, "allow_extra_args": True}


def _repo_root() -> Path:
    """Find the Ultidock repository root from an editable install or cwd."""
    candidates = [
        Path(__file__).resolve().parent.parent,
        Path.cwd(),
        Path.cwd().parent,
    ]
    for candidate in candidates:
        if (
            (candidate / "docking").is_dir()
            and (candidate / "benchmarks").is_dir()
            and (candidate / "molguard").is_dir()
        ):
            return candidate
    return Path.cwd()


def _require_dir(path: Path, label: str, *, topic: str = "quick-start") -> Path:
    if not path.is_dir():
        click.echo(f"[FAIL] Could not find {label}: {path}", err=True)
        click.echo(f"       {readme_hint(topic)}", err=True)
        raise SystemExit(1)
    return path


def _require_file(path: Path, label: str, *, topic: str = "commands") -> Path:
    if not path.is_file():
        click.echo(f"[FAIL] Could not find {label}: {path}", err=True)
        click.echo(f"       {readme_hint(topic)}", err=True)
        raise SystemExit(1)
    return path


def _docking_dir() -> Path:
    return _require_dir(_repo_root() / "docking", "docking directory", topic="quick-start")


def _benchmarks_dir() -> Path:
    return _require_dir(_repo_root() / "benchmarks", "benchmarks directory", topic="benchmarks")


def _examples_dir() -> Path:
    return _require_dir(_repo_root() / "examples", "examples directory", topic="examples")


def _run_python(
    script: Path,
    args: tuple[str, ...],
    *,
    cwd: Path | None = None,
    topic: str = "commands",
) -> None:
    script = _require_file(script, "script", topic=topic)
    result = subprocess.run([sys.executable, "-u", str(script), *args], cwd=cwd or _repo_root())
    if result.returncode:
        click.echo(f"       {readme_hint(topic)}", err=True)
    raise SystemExit(result.returncode)


def _find_tool(binary: str, extra_dirs: list[Path]) -> str | None:
    import shutil

    for directory in extra_dirs:
        candidate = directory / binary
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    return shutil.which(binary)


@click.group()
@click.version_option(version=__version__, prog_name="ultidock")
def cli() -> None:
    """Ultidock: run docking workflows, examples, and benchmarks."""


@cli.command("doctor")
def doctor_cmd() -> None:
    """Print Ultidock tool paths and environment diagnostics."""
    root = _repo_root()
    local_dirs: list[Path] = [
        root / "docking" / "AUTODOCK_GPU_DIR" / "autogrid",
        root / "docking" / "AUTODOCK_GPU_DIR" / "bin",
        root / "docking" / "VINA_DIR" / "bin",
        root / "docking" / ".toolshims",
    ]
    tools = [
        ("autogrid4", "AutoGrid", root / "docking" / "AUTODOCK_GPU_DIR" / "autogrid"),
        ("autodock4", "AutoDock 4", None),
        ("autodock_gpu_128wi", "AutoDock-GPU (128wi)", root / "docking" / "AUTODOCK_GPU_DIR"),
        ("autodock_gpu_64wi", "AutoDock-GPU (64wi)", root / "docking" / "AUTODOCK_GPU_DIR"),
        ("vina", "AutoDock Vina", None),
    ]

    click.echo("-- ultidock environment --------------------------------")
    click.echo(f"  {'ultidock':20s} {__version__}")
    click.echo(f"  {'python':20s} {sys.version.split()[0]}")
    click.echo(f"  {'repo root':20s} {root}")

    click.echo("-- external tools --------------------------------------")
    had_issue = False
    for binary, label, source_dir in tools:
        found = _find_tool(binary, local_dirs)
        if found:
            click.echo(f"  [OK]    {label:30s} {found}")
        elif source_dir is not None and source_dir.is_dir():
            had_issue = True
            click.echo(f"  [WARN]  {label:30s} not compiled  (source: {source_dir})", err=True)
        else:
            had_issue = True
            click.echo(f"  [FAIL]  {label:30s} not found", err=True)
    if had_issue:
        click.echo(f"       {readme_hint('troubleshooting')}", err=True)


@cli.command("setup", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def setup_cmd(extra_args: tuple[str, ...]) -> None:
    """Run docking/setup.py with all extra arguments forwarded."""
    _run_python(_docking_dir() / "setup.py", extra_args, cwd=_docking_dir(), topic="requirements")


@cli.command("run", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def run_cmd(extra_args: tuple[str, ...]) -> None:
    """Run the full Ultidock docking pipeline."""
    _run_python(_docking_dir() / "run.py", extra_args, cwd=_docking_dir(), topic="quick-start")


@cli.command("clean", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def clean_cmd(extra_args: tuple[str, ...]) -> None:
    """Reset compiled binaries and docking outputs."""
    _run_python(_docking_dir() / "clean.py", extra_args, cwd=_docking_dir(), topic="troubleshooting")


@cli.command("profile-receptors", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def profile_receptors_cmd(extra_args: tuple[str, ...]) -> None:
    """Generate optional per-receptor geometry sidecars."""
    _run_python(_docking_dir() / "profile_receptors.py", extra_args, cwd=_docking_dir(), topic="site-finder")


@cli.group()
def benchmark() -> None:
    """Benchmark download, site-recovery, docking, and plotting commands."""


@benchmark.command("download-dude", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_download_dude_cmd(extra_args: tuple[str, ...]) -> None:
    """Download DUD-E receptor, crystal ligand, active, and decoy files."""
    _run_python(_benchmarks_dir() / "download_dude.py", extra_args, topic="benchmarks")


@benchmark.command("cavity-recovery", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_cavity_recovery_cmd(extra_args: tuple[str, ...]) -> None:
    """Evaluate receptor-only site recovery against crystal-ligand centers."""
    _run_python(_benchmarks_dir() / "cavity_recovery_benchmark.py", extra_args, topic="benchmarks")


@benchmark.command("dude-docking", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_dude_docking_cmd(extra_args: tuple[str, ...]) -> None:
    """Run the per-target DUD-E docking benchmark."""
    _run_python(_benchmarks_dir() / "dude_docking_benchmark.py", extra_args, topic="benchmarks")


@benchmark.command("full-dude", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_full_dude_cmd(extra_args: tuple[str, ...]) -> None:
    """Run the autonomous DUD-E download/dock/score/plot benchmark."""
    _run_python(_benchmarks_dir() / "run_full_dude_benchmark.py", extra_args, topic="benchmarks")


@benchmark.command("plots", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_plots_cmd(extra_args: tuple[str, ...]) -> None:
    """Generate benchmark plots from existing summary files."""
    _run_python(_benchmarks_dir() / "dude_plots.py", extra_args, topic="benchmarks")


@cli.group(name="example")
def example_group() -> None:
    """List and run bundled example pipelines."""


def _example_roots() -> list[Path]:
    return sorted(path for path in _examples_dir().iterdir() if (path / "example-run.py").is_file())


@example_group.command("list")
def example_list_cmd() -> None:
    """List runnable bundled examples."""
    roots = _example_roots()
    if not roots:
        click.echo("No runnable examples found.")
        return
    for path in roots:
        click.echo(path.name)


@example_group.command("run", context_settings=FORWARD_CONTEXT)
@click.argument("name")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def example_run_cmd(name: str, extra_args: tuple[str, ...]) -> None:
    """Run a bundled example by directory name."""
    matches = [path for path in _example_roots() if path.name == name]
    if not matches:
        available = ", ".join(path.name for path in _example_roots()) or "none"
        click.echo(f"[FAIL] Unknown example: {name}", err=True)
        click.echo(f"Available examples: {available}", err=True)
        click.echo(f"       {readme_hint('examples')}", err=True)
        raise SystemExit(1)
    _run_python(matches[0] / "example-run.py", extra_args, cwd=matches[0], topic="examples")


cli.add_command(example_group, "examples")
