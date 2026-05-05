"""MolGuard command-line interface for deterministic molecular file handling."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import click

from cli.readme import readme_hint


@click.group()
@click.version_option(package_name="molguard")
def cli() -> None:
    """MolGuard: deterministic molecular I/O checks and repair tools."""


@cli.group()
def pdbqt() -> None:
    """PDBQT ligand/receptor validation and normalization."""


@pdbqt.command("check")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
def pdbqt_check_cmd(file: Path) -> None:
    """Lint a PDBQT file and report all issues."""
    from molguard.io.pdbqt import pdbqt_check

    report = pdbqt_check(file)

    for error in report.errors:
        click.echo(
            f"ERROR   line {error.line_no:5d}  [{error.column:10s}]  "
            f"{error.code}: {error.message}",
            err=True,
        )
    for warning in report.warnings:
        click.echo(
            f"WARNING line {warning.line_no:5d}  [{warning.column:10s}]  "
            f"{warning.code}: {warning.message}",
            err=True,
        )

    if report.errors:
        click.echo(
            f"\n[FAIL] {file.name}: {len(report.errors)} error(s), "
            f"{len(report.warnings)} warning(s)",
            err=True,
        )
        click.echo(f"       {readme_hint('troubleshooting')}", err=True)
        raise SystemExit(1)

    click.echo(f"[OK] {file.name}: OK  ({len(report.warnings)} warning(s))")


@pdbqt.command("normalize")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for the normalized PDBQT.",
)
def pdbqt_normalize_cmd(file: Path, output: Path) -> None:
    """Reformat numeric columns in a ligand PDBQT without changing torsion trees."""
    from molguard.io.pdbqt import LintError, pdbqt_normalize

    try:
        pdbqt_normalize(file, output)
    except LintError as exc:
        click.echo(f"[FAIL] {exc}", err=True)
        click.echo(f"       {readme_hint('troubleshooting')}", err=True)
        raise SystemExit(1) from exc

    click.echo(f"[OK] Normalized -> {output}")


@cli.group()
def receptor() -> None:
    """Deterministic receptor preparation and canonicalization."""


@receptor.command("canonicalize")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for the canonicalized receptor PDBQT.",
)
@click.option(
    "--timestamp",
    default=None,
    help="Fixed REMARK timestamp for byte-identical output. Defaults to current UTC time.",
)
def receptor_canonicalize_cmd(file: Path, output: Path, timestamp: str | None) -> None:
    """Canonicalize a receptor PDBQT deterministically."""
    from molguard.io.pdbqt import LintError, canonicalize_receptor

    try:
        digest = canonicalize_receptor(file, output, timestamp=timestamp)
    except LintError as exc:
        click.echo(f"[FAIL] {exc}", err=True)
        click.echo(f"       {readme_hint('receptor-inputs')}", err=True)
        raise SystemExit(1) from exc

    click.echo(f"[OK] Canonicalized receptor -> {output}")
    click.echo(f"  sha256: {digest}")


@receptor.command("prepare")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for the prepared receptor PDBQT.",
)
@click.option(
    "--prepare-command",
    help="External conversion command template. Use {input}, {output}, and optionally {seed}.",
)
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--timestamp", default="MOLGUARD", show_default=True)
def receptor_prepare_cmd(
    file: Path,
    output: Path,
    prepare_command: str | None,
    seed: int,
    timestamp: str,
) -> None:
    """Prepare PDB/PDBQT/MOL2 receptor input through the shared MolGuard path."""
    from molguard.io.receptor_prep import prepare_receptor_pdbqt

    try:
        digest = prepare_receptor_pdbqt(
            input_path=file,
            output_path=output,
            prepare_command=prepare_command,
            seed=seed,
            timestamp=timestamp,
        )
    except Exception as exc:
        click.echo(f"[FAIL] {exc}", err=True)
        click.echo(f"       {readme_hint('receptor-inputs')}", err=True)
        raise SystemExit(1) from exc

    click.echo(f"[OK] Prepared receptor -> {output}")
    click.echo(f"  sha256: {digest}")


@pdbqt.command("canonicalize-receptor", hidden=True)
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", required=True, type=click.Path(path_type=Path))
def legacy_canonicalize_receptor_cmd(file: Path, output: Path) -> None:
    """Compatibility alias for `molguard receptor canonicalize`."""
    receptor_canonicalize_cmd.callback(file=file, output=output, timestamp=None)


@pdbqt.command("prepare-receptor", hidden=True)
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", required=True, type=click.Path(path_type=Path))
@click.option("--prepare-command")
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--timestamp", default="MOLGUARD", show_default=True)
def legacy_prepare_receptor_cmd(
    file: Path,
    output: Path,
    prepare_command: str | None,
    seed: int,
    timestamp: str,
) -> None:
    """Compatibility alias for `molguard receptor prepare`."""
    receptor_prepare_cmd.callback(
        file=file,
        output=output,
        prepare_command=prepare_command,
        seed=seed,
        timestamp=timestamp,
    )


@cli.group()
def grids() -> None:
    """AutoGrid .fld and .map file validation."""


@grids.command("check")
@click.argument("fld", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--types",
    multiple=True,
    metavar="TYPE",
    help=(
        "Expected atom types such as --types C --types HD --types OA. "
        "If omitted, types are derived from the .fld file."
    ),
)
def grids_check_cmd(fld: Path, types: tuple[str, ...]) -> None:
    """Sanity-check a .fld file and all referenced .map files."""
    from molguard.grids.check import check_fld

    report = check_fld(fld, expected_types=list(types) if types else None)

    for error in report.errors:
        click.echo(f"ERROR   (fld)  {error}", err=True)
    for warning in report.warnings:
        click.echo(f"WARNING (fld)  {warning}", err=True)

    for result in report.map_results:
        status = "OK  " if result.ok else "FAIL"
        range_str = ""
        if result.sampled_min is not None and result.sampled_max is not None:
            range_str = f"  range=[{result.sampled_min:.3f}, {result.sampled_max:.3f}]"
        click.echo(f"  [{status}]  {result.atom_type:4s}  {result.map_path.name}{range_str}")
        for error in result.errors:
            click.echo(f"           ERROR   {error}", err=True)
        for warning in result.warnings:
            click.echo(f"           WARNING {warning}", err=True)

    if report.ok:
        click.echo(f"[OK] All maps OK ({len(report.map_results)} checked)")
        return

    n_fail = sum(1 for result in report.map_results if not result.ok)
    click.echo(
        f"[FAIL] {n_fail} map(s) failed + {len(report.errors)} fld-level error(s)",
        err=True,
    )
    click.echo(f"       {readme_hint('troubleshooting')}", err=True)
    raise SystemExit(1)


@cli.command("doctor")
def doctor_cmd() -> None:
    """Print MolGuard version and optional receptor-conversion backends."""
    from molguard import __version__

    click.echo("-- molguard environment --------------------------------")
    click.echo(f"  {'molguard':24s} {__version__}")
    click.echo(f"  {'python':24s} {sys.version.split()[0]}")

    tools = [
        ("mk_prepare_receptor.py", "Meeko receptor prep"),
        ("mk_prepare_receptor", "Meeko receptor prep"),
        ("obabel", "Open Babel"),
    ]
    click.echo("-- optional conversion tools ---------------------------")
    had_missing = False
    for binary, label in tools:
        found = shutil.which(binary)
        if found:
            click.echo(f"  [OK]    {label:24s} {found}")
        else:
            had_missing = True
            click.echo(f"  [WARN]  {label:24s} not found", err=True)
    if had_missing:
        click.echo(f"       {readme_hint('receptor-inputs')}", err=True)
