"""Analysis exports must keep scores associated with their docking sites."""

from __future__ import annotations

import csv
import runpy
import sqlite3
import sys
from pathlib import Path
from types import ModuleType

import pytest


ANALYSIS_SCRIPT = Path(__file__).resolve().parents[2] / "docking/analyse_docking_results.py"
EXPORT_COLUMNS = [
    "ligand_name",
    "affinity",
    "rmsd_lb",
    "rmsd_ub",
    "docking_file",
    "binding_site",
    "model",
    "zinc_id",
]


def _export(monkeypatch, tmp_path, rows, *, legacy=False):
    database = tmp_path / "results.db"
    output = tmp_path / "results.csv"
    site_column = "" if legacy else ", binding_site TEXT"
    with sqlite3.connect(database) as connection:
        connection.execute(
            "CREATE TABLE docking_results ("
            'ligand_name TEXT, "binding_affinity (kcal/mol)" REAL, '
            '"rmsd_lb (Å)" REAL, "rmsd_ub (Å)" REAL, docking_file TEXT'
            f"{site_column})"
        )
        placeholders = ", ".join("?" for _ in range(5 if legacy else 6))
        connection.executemany(f"INSERT INTO docking_results VALUES ({placeholders})", rows)

    config = ModuleType("config")
    config.DB_PATH = str(database)
    monkeypatch.setitem(sys.modules, "config", config)
    monkeypatch.setattr(sys, "argv", [str(ANALYSIS_SCRIPT), "--out", str(output)])
    runpy.run_path(str(ANALYSIS_SCRIPT), run_name="__main__")
    with output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


@pytest.mark.parametrize("rmsds", [(0.0, 0.0), (None, None)], ids=["vina", "gpu"])
def test_export_preserves_sites_for_identical_pose_names(monkeypatch, tmp_path, rmsds):
    name = "escitalopram-e-Model1-5i6x_edited"
    ligand = "/inputs/escitalopram-e.pdbqt"
    rows = [
        (name, -8.5, *rmsds, ligand, "S2"),
        (name, -9.0, *rmsds, ligand, "1"),
        (name, -8.0, *rmsds, ligand, None),
        (name, -6.0, *rmsds, ligand, "S3"),
        (name, -10.0, 6.0, 0.0, ligand, "S4"),
        (name, -10.0, 0.0, 11.0, ligand, "S5"),
        ("escitalopram-e-Model0-5i6x_edited", -10.0, *rmsds, ligand, "S6"),
    ]

    columns, exported = _export(monkeypatch, tmp_path, rows)

    assert columns == EXPORT_COLUMNS
    assert [(row["binding_site"], float(row["affinity"])) for row in exported] == [
        ("1", -9.0),
        ("S2", -8.5),
        ("", -8.0),
    ]
    assert all(row["ligand_name"] == name for row in exported)
    assert all(row["docking_file"] == ligand for row in exported)


def test_empty_export_still_has_site_column(monkeypatch, tmp_path):
    columns, exported = _export(monkeypatch, tmp_path, [])

    assert columns == EXPORT_COLUMNS
    assert exported == []


def test_legacy_database_exports_unknown_site_without_guessing(monkeypatch, tmp_path):
    rows = [("ligand-Model1-receptor", -9.0, None, None, "/inputs/ligand.pdbqt")]

    columns, exported = _export(monkeypatch, tmp_path, rows, legacy=True)

    assert columns == EXPORT_COLUMNS
    assert exported[0]["binding_site"] == ""
