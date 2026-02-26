"""
ultidock.cli
~~~~~~~~~~~~
Command-line interface for Ultidock I/O hardening tools.

Commands:
    ultidock pdbqt check                <file.pdbqt>
    ultidock pdbqt normalize            <file.pdbqt> -o <out.pdbqt>
    ultidock pdbqt canonicalize-receptor <file.pdbqt> -o <out.pdbqt>
    ultidock grids check                <maps.fld>   [--types C HD OA ...]
    ultidock doctor

All commands print errors to stderr and use exit code 1 on failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

# ── Root group ────────────────────────────────────────────────────────────────

@click.group()
@click.version_option(package_name="molguard")
def cli() -> None:
    """Ultidock — GPU-accelerated docking pipeline tools."""


# ── pdbqt sub-group ───────────────────────────────────────────────────────────

@cli.group()
def pdbqt() -> None:
    """PDBQT file validation and normalization."""


@pdbqt.command("check")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
def pdbqt_check_cmd(file: Path) -> None:
    """Lint a PDBQT file and report all issues.

    Exits with code 0 on success (warnings allowed), code 1 if errors found.
    """
    from molguard.io.pdbqt import pdbqt_check

    report = pdbqt_check(file)

    for e in report.errors:
        click.echo(
            f"ERROR   line {e.line_no:5d}  [{e.column:10s}]  {e.code}: {e.message}",
            err=True,
        )
    for w in report.warnings:
        click.echo(
            f"WARNING line {w.line_no:5d}  [{w.column:10s}]  {w.code}: {w.message}",
            err=True,
        )

    if report.errors:
        click.echo(
            f"\n✗ {file.name}: {len(report.errors)} error(s), "
            f"{len(report.warnings)} warning(s)",
            err=True,
        )
        sys.exit(1)

    click.echo(
        f"✓ {file.name}: OK  ({len(report.warnings)} warning(s))"
    )


@pdbqt.command("normalize")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for the normalized PDBQT.",
)
def pdbqt_normalize_cmd(file: Path, output: Path) -> None:
    """Reformat numeric columns in a ligand PDBQT in-place.

    Torsion tree structure (ROOT/BRANCH/ENDBRANCH/TORSDOF) is preserved
    verbatim. Serial numbers and atom order are not changed.
    """
    from molguard.io.pdbqt import LintError, pdbqt_normalize

    try:
        pdbqt_normalize(file, output)
    except LintError as exc:
        click.echo(f"✗ {exc}", err=True)
        sys.exit(1)

    click.echo(f"✓ Normalized → {output}")


@pdbqt.command("canonicalize-receptor")
@click.argument("file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output path for the canonicalized receptor PDBQT.",
)
def canonicalize_cmd(file: Path, output: Path) -> None:
    """Deterministically canonicalize a receptor PDBQT.

    Atoms are sorted by chain → residue → atom name, serials renumbered,
    and all numeric fields reformatted through the fixed-width formatter.
    Running twice on the same input produces byte-identical output.
    """
    from molguard.io.pdbqt import LintError, canonicalize_receptor

    try:
        digest = canonicalize_receptor(file, output)
    except LintError as exc:
        click.echo(f"✗ {exc}", err=True)
        sys.exit(1)

    click.echo(f"✓ Canonicalized → {output}")
    click.echo(f"  sha256: {digest}")


# ── grids sub-group ───────────────────────────────────────────────────────────

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
        "Expected atom types separated by spaces (e.g. --types C HD OA N e d). "
        "If omitted, types are derived from the .fld file itself and no "
        "cross-check against expected types is performed."
    ),
)
def grids_check_cmd(fld: Path, types: tuple[str, ...]) -> None:
    """Sanity-check a .fld file and all referenced .map files.

    Detects all-zero maps (AutoGrid read failure), NaN/Inf, absurd energy
    ranges, missing .map files, and atom-type count mismatches.
    """
    from molguard.grids.check import check_fld

    report = check_fld(fld, expected_types=list(types) if types else None)

    # fld-level errors
    for e in report.errors:
        click.echo(f"ERROR   (fld)  {e}", err=True)
    for w in report.warnings:
        click.echo(f"WARNING (fld)  {w}", err=True)

    # per-map results
    for r in report.map_results:
        status = "OK  " if r.ok else "FAIL"
        range_str = ""
        if r.sampled_min is not None and r.sampled_max is not None:
            range_str = f"  range=[{r.sampled_min:.3f}, {r.sampled_max:.3f}]"
        click.echo(f"  [{status}]  {r.atom_type:4s}  {r.map_path.name}{range_str}")
        for e in r.errors:
            click.echo(f"           ERROR   {e}", err=True)
        for w in r.warnings:
            click.echo(f"           WARNING {w}", err=True)

    if report.ok:
        click.echo(f"✓ All maps OK ({len(report.map_results)} checked)")
        sys.exit(0)
    else:
        n_fail = sum(1 for r in report.map_results if not r.ok)
        click.echo(
            f"✗ {n_fail} map(s) failed + {len(report.errors)} fld-level error(s)",
            err=True,
        )
        sys.exit(1)


# ── doctor ────────────────────────────────────────────────────────────────────

@cli.command("doctor")
def doctor_cmd() -> None:
    """Print tool paths, versions, and environment diagnostics."""
    import shutil

    from molguard import __version__

    click.echo("── ultidock environment ──────────────────────────────────")
    click.echo(f"  {'ultidock':20s} {__version__}")
    click.echo(f"  {'python':20s} {sys.version.split()[0]}")

    tools = [
        ("autogrid4",           "AutoGrid"),
        ("autodock4",           "AutoDock 4"),
        ("autodock_gpu_128wi",  "AutoDock-GPU (128wi)"),
        ("vina",                "AutoDock Vina"),
    ]
    click.echo("── external tools ────────────────────────────────────────")
    for binary, label in tools:
        p = shutil.which(binary)
        if p:
            click.echo(f"  ✓  {label:30s} {p}")
        else:
            click.echo(f"  ✗  {label:30s} not found in PATH", err=True)
