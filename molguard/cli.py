"""
molguard.cli
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
            f"\n[FAIL] {file.name}: {len(report.errors)} error(s), "
            f"{len(report.warnings)} warning(s)",
            err=True,
        )
        sys.exit(1)

    click.echo(f"[OK] {file.name}: OK  ({len(report.warnings)} warning(s))")


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
        click.echo(f"[FAIL] {exc}", err=True)
        sys.exit(1)

    click.echo(f"[OK] Normalized -> {output}")


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

    Atoms are sorted by chain -> residue -> atom name, serials renumbered,
    and all numeric fields reformatted through the fixed-width formatter.
    Running twice on the same input produces byte-identical output.
    """
    from molguard.io.pdbqt import LintError, canonicalize_receptor

    try:
        digest = canonicalize_receptor(file, output)
    except LintError as exc:
        click.echo(f"[FAIL] {exc}", err=True)
        sys.exit(1)

    click.echo(f"[OK] Canonicalized -> {output}")
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
        "Expected atom types (e.g. --types C HD OA N e d). "
        "If omitted, types are derived from the .fld file and no cross-check is performed."
    ),
)
def grids_check_cmd(fld: Path, types: tuple[str, ...]) -> None:
    """Sanity-check a .fld file and all referenced .map files.

    Detects all-zero maps (AutoGrid read failure), NaN/Inf, absurd energy
    ranges, missing .map files, and atom-type count mismatches.
    """
    from molguard.grids.check import check_fld

    report = check_fld(fld, expected_types=list(types) if types else None)

    for e in report.errors:
        click.echo(f"ERROR   (fld)  {e}", err=True)
    for w in report.warnings:
        click.echo(f"WARNING (fld)  {w}", err=True)

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
        click.echo(f"[OK] All maps OK ({len(report.map_results)} checked)")
        sys.exit(0)
    else:
        n_fail = sum(1 for r in report.map_results if not r.ok)
        click.echo(
            f"[FAIL] {n_fail} map(s) failed + {len(report.errors)} fld-level error(s)",
            err=True,
        )
        sys.exit(1)


# ── doctor ────────────────────────────────────────────────────────────────────

def _repo_root() -> Path:
    """
    Walk upward from this file to find the repo root, identified by the
    presence of both 'docking/' and 'molguard/' directories.
    Falls back to cwd if the heuristic fails.
    """
    for candidate in [
        Path(__file__).resolve().parent.parent,  # <repo>/molguard/../  = <repo>
        Path.cwd(),
        Path.cwd().parent,
    ]:
        if (candidate / "docking").is_dir() and (candidate / "molguard").is_dir():
            return candidate
    return Path.cwd()


def _find_tool(binary: str, extra_dirs: list[Path]) -> str | None:
    """
    Check extra_dirs first (repo-local), then fall back to system PATH.
    Returns the resolved absolute path string, or None if not found.
    """
    import os
    import shutil

    for d in extra_dirs:
        candidate = d / binary
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    return shutil.which(binary)


@cli.command("doctor")
def doctor_cmd() -> None:
    """Print tool paths, versions, and environment diagnostics."""
    from molguard import __version__

    root = _repo_root()

    # Directories to probe (in order) before falling back to system PATH.
    # Mirrors the layout produced by docking/setup.py / autodock-gpu-compiler.sh.
    local_dirs: list[Path] = [
        root / "docking" / "AUTODOCK_GPU_DIR" / "autogrid",  # autogrid4
        root / "docking" / "AUTODOCK_GPU_DIR" / "bin",       # autodock_gpu_*
        root / "docking" / "VINA_DIR" / "bin",               # vina, vina_split
        root / "docking" / ".toolshims",                      # any shim wrappers
    ]

    tools = [
        # (binary,            label,                   source_dir_hint)
        ("autogrid4",
         "AutoGrid",
         root / "docking" / "AUTODOCK_GPU_DIR" / "autogrid"),
        ("autodock4",
         "AutoDock 4",
         None),
        ("autodock_gpu_128wi",
         "AutoDock-GPU (128wi)",
         root / "docking" / "AUTODOCK_GPU_DIR"),
        ("autodock_gpu_64wi",
         "AutoDock-GPU (64wi)",
         root / "docking" / "AUTODOCK_GPU_DIR"),
        ("vina",
         "AutoDock Vina",
         None),
    ]

    click.echo("── ultidock environment ──────────────────────────────────")
    click.echo(f"  {'ultidock':20s} {__version__}")
    click.echo(f"  {'python':20s} {sys.version.split()[0]}")
    click.echo(f"  {'repo root':20s} {root}")

    click.echo("── external tools ────────────────────────────────────────")
    for binary, label, src_dir in tools:
        p = _find_tool(binary, local_dirs)
        if p:
            click.echo(f"  [OK]  {label:30s} {p}")
        elif src_dir is not None and src_dir.is_dir():
            click.echo(f"  [WARN]  {label:30s} not compiled  (source: {src_dir})", err=True)
        else:
            click.echo(f"  [FAIL]  {label:30s} not found", err=True)
