"""Ultidock command-line interface for workflow orchestration."""

from __future__ import annotations

import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import click

from cli.report import generate_report
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


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _default_run_dir(mode: str) -> Path:
    return (_repo_root() / "runs" / f"{mode}_{_timestamp()}").resolve()


def _split_center(center: str) -> tuple[float, float, float]:
    parts = [part.strip() for part in center.replace(";", ",").split(",")]
    if len(parts) != 3:
        raise click.BadParameter("center must be formatted as x,y,z")
    try:
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    except ValueError as exc:
        raise click.BadParameter("center must contain three numbers") from exc


def _write_run_config(path: Path, rows: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for key, value in rows.items():
        if isinstance(value, (list, tuple)):
            rendered = "[" + ", ".join(str(item) for item in value) + "]"
        else:
            rendered = str(value)
        lines.append(f"{key}: {rendered}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_known_site_tsv(path: Path, center: tuple[float, float, float], box_size: float) -> None:
    spacing = 0.375
    npts = max(1, int(round(float(box_size) / spacing)))
    if npts % 2 == 0:
        npts += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tr_peak\tF\t"
        "raw_F\tfamily\tportfolio_role\tselection_score\tcenter_closeness\n"
        "# meta policy=known_site site_count=1 ranking_score=manual_box\n"
        f"known_site\tS1\t{center[0]:.3f}\t{center[1]:.3f}\t{center[2]:.3f}\t"
        f"{npts}\t{npts}\t{npts}\t{spacing:.3f}\t{float(box_size) / 4.0:.2f}\t"
        "1.000\t1.000\tknown_site\tknown_site\t1.0\t1.0\n",
        encoding="utf-8",
    )


def _pipeline_dirs(run_dir: Path) -> dict[str, Path]:
    return {
        "docking": run_dir / "docking",
        "analysis": run_dir / "analysis",
        "results": run_dir / "results",
    }


def _run_pipeline_mode(
    *,
    public_mode: str,
    grid_mode: str,
    output_dir: Path | None,
    extra_args: tuple[str, ...],
    dry_run: bool,
    report: bool,
    centers_tsv: Path | None = None,
    config_extra: dict[str, object] | None = None,
) -> None:
    run_dir = (output_dir.resolve() if output_dir else _default_run_dir(public_mode))
    dirs = _pipeline_dirs(run_dir)
    for path in [run_dir, *dirs.values()]:
        path.mkdir(parents=True, exist_ok=True)
    if centers_tsv is None:
        centers_tsv = run_dir / "sites.tsv"

    command = [
        sys.executable,
        "-u",
        str(_docking_dir() / "run.py"),
        "--grid-mode",
        grid_mode,
        "--centers-tsv",
        str(centers_tsv),
        "--docking-dir",
        str(dirs["docking"]),
        "--analysis-dir",
        str(dirs["analysis"]),
        "--results-dir",
        str(dirs["results"]),
        *extra_args,
    ]
    config = {
        "workflow": public_mode,
        "site_method": "cav-emps" if public_mode == "cavity" else public_mode,
        "grid_mode": grid_mode,
        "run_dir": str(run_dir),
        "centers_tsv": str(centers_tsv),
        "command": " ".join(command),
    }
    if config_extra:
        config.update(config_extra)
    _write_run_config(run_dir / "run_config.yaml", config)

    if dry_run:
        click.echo("Dry run command:")
        click.echo(" ".join(command))
        click.echo(f"Run config: {run_dir / 'run_config.yaml'}")
        return

    result = subprocess.run(command, cwd=_docking_dir())
    if report:
        os.environ["ULTIDOCK_COMMAND"] = " ".join(command)
        generate_report(run_dir)
        click.echo(f"Report: {run_dir / 'report.html'}")
    if result.returncode:
        click.echo(f"       {readme_hint('quick-start')}", err=True)
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


@cli.command("known-site", context_settings=FORWARD_CONTEXT)
@click.option("--center", required=True, help="Known binding-site center as x,y,z.")
@click.option("--box-size", default=35.0, show_default=True, help="Manual docking box side in A.")
@click.option("--output-dir", type=click.Path(path_type=Path), help="Run directory.")
@click.option("--dry-run", is_flag=True, help="Write config and print the pipeline command only.")
@click.option("--report/--no-report", default=True, show_default=True, help="Generate report files.")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def known_site_cmd(
    center: str,
    box_size: float,
    output_dir: Path | None,
    dry_run: bool,
    report: bool,
    extra_args: tuple[str, ...],
) -> None:
    """Run docking with an expert/manual known-site box."""
    run_dir = output_dir.resolve() if output_dir else _default_run_dir("known_site")
    centers_tsv = run_dir / "sites.tsv"
    parsed_center = _split_center(center)
    _write_known_site_tsv(centers_tsv, parsed_center, box_size)
    _run_pipeline_mode(
        public_mode="known-site",
        grid_mode="centers",
        output_dir=run_dir,
        extra_args=extra_args,
        dry_run=dry_run,
        report=report,
        centers_tsv=centers_tsv,
        config_extra={"known_center": parsed_center, "box_size_a": box_size},
    )


@cli.command("cavity", context_settings=FORWARD_CONTEXT)
@click.option("--autosites", default=6, show_default=True, help="Number of CaV-EMPS sites.")
@click.option("--output-dir", type=click.Path(path_type=Path), help="Run directory.")
@click.option("--dry-run", is_flag=True, help="Write config and print the pipeline command only.")
@click.option("--report/--no-report", default=True, show_default=True, help="Generate report files.")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def cavity_cmd(
    autosites: int,
    output_dir: Path | None,
    dry_run: bool,
    report: bool,
    extra_args: tuple[str, ...],
) -> None:
    """Run CaV-EMPS automatic site proposal, then docking."""
    _run_pipeline_mode(
        public_mode="cavity",
        grid_mode="centers",
        output_dir=output_dir,
        extra_args=("--autosites", str(autosites), *extra_args),
        dry_run=dry_run,
        report=report,
        config_extra={"autosites": autosites},
    )


@cli.command("blind", context_settings=FORWARD_CONTEXT)
@click.option("--grid-cap", default=150.0, show_default=True, help="Blind box cap in A.")
@click.option("--output-dir", type=click.Path(path_type=Path), help="Run directory.")
@click.option("--dry-run", is_flag=True, help="Write config and print the pipeline command only.")
@click.option("--report/--no-report", default=True, show_default=True, help="Generate report files.")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def blind_cmd(
    grid_cap: float,
    output_dir: Path | None,
    dry_run: bool,
    report: bool,
    extra_args: tuple[str, ...],
) -> None:
    """Run naive whole-receptor/blind docking."""
    _run_pipeline_mode(
        public_mode="blind",
        grid_mode="blind",
        output_dir=output_dir,
        extra_args=("--grid-cap", str(grid_cap), *extra_args),
        dry_run=dry_run,
        report=report,
        config_extra={"grid_cap_a": grid_cap},
    )


@cli.command("report")
@click.argument("run_dir", type=click.Path(path_type=Path))
def report_cmd(run_dir: Path) -> None:
    """Generate Markdown/HTML reports and visualization files for a run directory."""
    outputs = generate_report(run_dir)
    click.echo(f"Markdown: {outputs['report_md']}")
    click.echo(f"HTML:     {outputs['report_html']}")


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


@benchmark.command("baseline-p2rank", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_baseline_p2rank_cmd(extra_args: tuple[str, ...]) -> None:
    """Evaluate P2Rank site recovery against crystal-ligand centers."""
    _run_python(
        _benchmarks_dir() / "baseline" / "p2rank" / "run_p2rank_baseline.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("baseline-fpocket", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_baseline_fpocket_cmd(extra_args: tuple[str, ...]) -> None:
    """Evaluate fpocket site recovery against crystal-ligand centers."""
    _run_python(
        _benchmarks_dir() / "baseline" / "fpocket" / "run_fpocket_baseline.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("site-evaluate", context_settings=FORWARD_CONTEXT, hidden=True)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_site_evaluate_cmd(extra_args: tuple[str, ...]) -> None:
    """Evaluate normalized site-prediction outputs with DCC Top-n metrics."""
    _run_python(
        _benchmarks_dir() / "site_prediction" / "evaluate_predictions.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("site-prediction", context_settings=FORWARD_CONTEXT)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_site_prediction_cmd(extra_args: tuple[str, ...]) -> None:
    """Download, normalize, run, evaluate, and report site-prediction benchmarks."""
    _run_python(
        _benchmarks_dir() / "site_prediction" / "run_site_prediction_benchmark.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("site-import", context_settings=FORWARD_CONTEXT, hidden=True)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_site_import_cmd(extra_args: tuple[str, ...]) -> None:
    """Advanced: download/normalize COACH420/HOLO4K datasets."""
    _run_python(
        _benchmarks_dir() / "site_prediction" / "prepare_site_datasets.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("site-run", context_settings=FORWARD_CONTEXT, hidden=True)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_site_run_cmd(extra_args: tuple[str, ...]) -> None:
    """Advanced: generate cav-emps, fpocket, or p2rank predictions."""
    _run_python(
        _benchmarks_dir() / "site_prediction" / "predict_sites.py",
        extra_args,
        topic="benchmarks",
    )


@benchmark.command("site-report", context_settings=FORWARD_CONTEXT, hidden=True)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def benchmark_site_report_cmd(extra_args: tuple[str, ...]) -> None:
    """Generate a Markdown/HTML site-benchmark summary page."""
    _run_python(
        _benchmarks_dir() / "site_prediction" / "report.py",
        extra_args,
        topic="benchmarks",
    )


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
    roots = sorted(roots, key=lambda path: (path.name != "quickstart", path.name))
    for path in roots:
        suffix = " (recommended first run)" if path.name == "quickstart" else ""
        click.echo(f"{path.name}{suffix}")


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
