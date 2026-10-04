#!/usr/bin/env python3
"""Helper utilities for running Ultidock example pipelines."""

from __future__ import annotations

import shutil
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_DIR = REPO_ROOT / "docking"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def add_pipeline_options(parser):
    from docking.run import add_md_arguments
    add_md_arguments(parser)
    parser.add_argument("--output-dir", type=Path, help="New, isolated example workspace")
    parser.add_argument("--dry-run", action="store_true", help="Stage inputs and print the command without running tools")


def pipeline_options(args):
    forwarded = []
    for name in ("md_config", "md_through", "md_work_dir", "md_gmx", "md_acpype", "md_obabel"):
        value = getattr(args, name)
        if value is not None:
            if name in {"md_config", "md_work_dir"} or (name.startswith("md_") and "/" in str(value)):
                value = Path(value).expanduser().resolve()
            forwarded.extend(["--" + name.replace("_", "-"), str(value)])
    return {**({"extra_args": forwarded} if forwarded else {}),
            **({"dry_run": True} if args.dry_run else {})}


def _ensure_workspace(example_root: Path, output_dir: Path | None = None) -> dict[str, Path]:
    run_name = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    workspace = output_dir.expanduser().resolve() if output_dir else example_root / "workspace" / run_name
    workspace.mkdir(parents=True, exist_ok=False)
    paths = {
        "workspace": workspace,
        "ligands": workspace / "LIGANDS_DIR",
        "macro": workspace / "MACRO_MOL_DIR",
        "docking": workspace / "DOCKING_DIR",
        "analysis": workspace / "ANALYSIS_DIR",
        "results": workspace / "RESULTS_DIR",
        "vina": DOCKING_DIR / "VINA_DIR",
        "autodock": DOCKING_DIR / "AUTODOCK_GPU_DIR",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def stage_inputs(example_root: Path, receptor: Path, ligands: Iterable[Path], *, output_dir=None) -> dict[str, Path]:
    """Prepare an isolated workspace for an example run."""

    paths = _ensure_workspace(example_root, output_dir)
    print(f"Example workspace: {paths['workspace']}")
    receptor_target = paths["macro"] / receptor.name
    shutil.copy2(receptor, receptor_target)
    paths["receptor"] = receptor_target

    for ligand in ligands:
        shutil.copy2(ligand, paths["ligands"] / ligand.name)

    # The binaries are shared; the inputs and results stay in this run's workspace.
    (paths["autodock"] / "bin").mkdir(parents=True, exist_ok=True)
    (paths["autodock"] / "autogrid").mkdir(parents=True, exist_ok=True)
    (paths["vina"] / "bin").mkdir(parents=True, exist_ok=True)

    return paths


def run_pipeline(
    paths: Mapping[str, Path], *, mode: str = "cpu", site_method: str = "cav-emps",
    extra_args=(), dry_run=False,
) -> None:
    """Execute the full Ultidock pipeline using the staged workspace."""

    if site_method not in {"cav-emps", "p2rank", "fpocket"}:
        raise ValueError(f"Unsupported site method: {site_method}")

    run_cmd = [
        sys.executable,
        "-u",
        "run.py",
        "--mode",
        mode,
        "--skip-wget",
        "--LIGANDS_DIR",
        str(paths["ligands"].resolve()),
        "--DOCKING_DIR",
        str(paths["docking"].resolve()),
        "--ANALYSIS_DIR",
        str(paths["analysis"].resolve()),
        "--VINA_DIR",
        str(paths["vina"].resolve()),
        "--AUTODOCK_GPU_DIR",
        str(paths["autodock"].resolve()),
        "--MACRO_MOL_DIR",
        str(paths["macro"].resolve()),
        "--RESULTS_DIR",
        str(paths["results"].resolve()),
    ]

    centers_tsv = paths["results"] / f"{site_method}-sites.tsv"
    if site_method in {"p2rank", "fpocket"} and not dry_run:
        if str(REPO_ROOT) not in sys.path:
            sys.path.insert(0, str(REPO_ROOT))
        from docking.pocket_boxes import create_pocket_boxes

        create_pocket_boxes(
            method=site_method,
            receptor_pdbqt=paths["receptor"],
            output_tsv=centers_tsv,
            work_dir=paths["results"] / f"{site_method}-pockets",
        )

    run_cmd.extend(["--grid-mode", "centers", "--centers-tsv", str(centers_tsv.resolve())])
    run_cmd.extend(extra_args)
    if dry_run:
        print(f"Site method: {site_method}")
        print(f"Pipeline command (from {DOCKING_DIR}): {shlex.join(run_cmd)}")
        return

    subprocess.run(run_cmd, check=True, cwd=DOCKING_DIR)
