#!/usr/bin/env python3
"""Helper utilities for running Ultidock example pipelines."""

from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_DIR = REPO_ROOT / "docking"


def _ensure_workspace(example_root: Path) -> dict[str, Path]:
    run_name = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    workspace = example_root / "workspace" / run_name
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


def stage_inputs(example_root: Path, receptor: Path, ligands: Iterable[Path]) -> dict[str, Path]:
    """Prepare an isolated workspace for an example run."""

    paths = _ensure_workspace(example_root)
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
    paths: Mapping[str, Path], *, mode: str = "cpu", site_method: str = "cav-emps"
) -> None:
    """Execute the full Ultidock pipeline using the staged workspace."""

    if site_method not in {"cav-emps", "p2rank", "fpocket"}:
        raise ValueError(f"Unsupported site method: {site_method}")

    run_cmd = [
        sys.executable,
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
    if site_method in {"p2rank", "fpocket"}:
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

    subprocess.run(run_cmd, check=True, cwd=DOCKING_DIR)
