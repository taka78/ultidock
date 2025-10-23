#!/usr/bin/env python3
"""Helper utilities for running Ultidock example pipelines."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCKING_DIR = REPO_ROOT / "docking" #This can be confusing but do i look like i care?


def _ensure_workspace(example_root: Path) -> Mapping[str, Path]:
    workspace = example_root / "workspace"
    paths = {
        "workspace": DOCKING_DIR,
        "ligands": DOCKING_DIR / "LIGANDS_DIR",
        "macro": DOCKING_DIR / "MACRO_MOL_DIR",
        "docking": DOCKING_DIR / "DOCKING_DIR",
        "analysis": DOCKING_DIR / "ANALYSIS_DIR",
        "results": DOCKING_DIR / "RESULTS_DIR",
        "vina": DOCKING_DIR / "VINA_DIR",
        "autodock": DOCKING_DIR / "AUTODOCK_GPU_DIR",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def stage_inputs(example_root: Path, receptor: Path, ligands: Iterable[Path]) -> Mapping[str, Path]:
    """Prepare an isolated workspace for an example run."""

    paths = _ensure_workspace(example_root)
    receptor_target = paths["macro"] / receptor.name
    shutil.copy2(receptor, receptor_target)

    for ligand in ligands:
        shutil.copy2(ligand, paths["ligands"] / ligand.name)

    # Provide placeholder bin directories the setup script expects.
    (paths["autodock"] / "bin").mkdir(parents=True, exist_ok=True)
    (paths["autodock"] / "autogrid").mkdir(parents=True, exist_ok=True)
    (paths["vina"] / "bin").mkdir(parents=True, exist_ok=True)

    return paths


def run_pipeline(paths: Mapping[str, Path], *, mode: str = "cpu") -> None:
    """Execute the full Ultidock pipeline using the staged workspace."""

    run_cmd = [
        "python3",
        "run.py",
        "--mode",
        mode,
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

    subprocess.run(run_cmd, check=True, cwd=DOCKING_DIR)
