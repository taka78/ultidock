#!/usr/bin/env python3
"""Generate the official Ultidock quickstart run artifacts."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from cli.report import generate_report  # noqa: E402
from molguard import __version__  # noqa: E402

FIXTURE_DIR = Path(__file__).resolve().parent / "data"


def _stage_fixtures(output_dir: Path) -> tuple[Path, Path]:
    fixture_output = output_dir / "input"
    fixture_output.mkdir(parents=True, exist_ok=True)
    receptor = fixture_output / "receptor.pdb"
    ligand = fixture_output / "reference_ligand.mol2"
    shutil.copy2(FIXTURE_DIR / "receptor.pdb", receptor)
    shutil.copy2(FIXTURE_DIR / "reference_ligand.mol2", ligand)
    return receptor, ligand


def _write_config(path: Path, *, receptor: Path, reference_ligand: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "workflow: quickstart",
                "site_method: cav-emps",
                "execution_mode: bundled_quickstart_fixture",
                f"version: {__version__}",
                f"command: {shlex.join([sys.executable, *sys.argv])}",
                "description: lightweight first-run artifact generation",
                f"receptor_input: {receptor}",
                f"reference_ligand: {reference_ligand}",
                "reference_ligand_used_for_site_generation: false",
                f"generated_at_utc: {datetime.now(timezone.utc).isoformat()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_sites(path: Path) -> None:
    path.write_text(
        "# receptor\tsite_id\tcx\tcy\tcz\tnx\tny\tnz\tspacing\tr_peak\tF\t"
        "raw_F\tfamily\tportfolio_role\tselection_score\tcenter_closeness\n"
        "# meta policy=quickstart site_count=3 ranking_score=cav-emps_example\n"
        "quickstart_receptor\tS1\t12.500\t8.250\t-3.000\t95\t95\t95\t0.375\t4.20\t"
        "1.080\t0.240\tsurface\tsurface_focus\t0.91\t0.82\n"
        "quickstart_receptor\tS2\t-4.250\t10.750\t6.500\t95\t95\t95\t0.375\t5.10\t"
        "0.990\t0.000\tcore_reserve\tcore_reserve\t\t\n"
        "quickstart_receptor\tS3\t21.000\t-2.500\t4.750\t95\t95\t95\t0.375\t3.40\t"
        "0.720\t0.310\thybrid\thybrid_primary\t0.65\t0.40\n",
        encoding="utf-8",
    )


def _write_predictions(path: Path) -> None:
    rows = [
        ["quickstart", "cav-emps", 1, "S1", 12.5, 8.25, -3.0, 1.08, "sites.tsv"],
        ["quickstart", "cav-emps", 2, "S2", -4.25, 10.75, 6.5, 0.99, "sites.tsv"],
        ["quickstart", "cav-emps", 3, "S3", 21.0, -2.5, 4.75, 0.72, "sites.tsv"],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["target_id", "method", "rank", "site_id", "center_x", "center_y", "center_z", "score", "source"]
        )
        writer.writerows(rows)


def _write_top_hits(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand_id", "site_id", "score", "pose_path", "note"])
        writer.writerow(["LIG001", "S1", "-8.2", "poses/LIG001_S1.pdbqt", "example row"])
        writer.writerow(["LIG002", "S1", "-7.6", "poses/LIG002_S1.pdbqt", "example row"])


def _write_sqlite(path: Path) -> None:
    with sqlite3.connect(path) as conn:
        conn.execute("create table if not exists metadata (key text primary key, value text)")
        conn.execute(
            "insert or replace into metadata values (?, ?)",
            ("workflow", "quickstart bundled fixture"),
        )
        conn.commit()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create the official Ultidock quickstart artifact directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        default="workspace/quickstart_run",
        help="Directory where quickstart artifacts will be written.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned outputs and exit.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir).resolve()
    if args.dry_run:
        print(f"Quickstart output dir: {output_dir}")
        print(
            "Would write: input/receptor.pdb, input/reference_ligand.mol2, "
            "run_config.yaml, sites.tsv, predictions.tsv, top_hits.csv, "
            "results.sqlite, report.md, report.html"
        )
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "poses").mkdir(exist_ok=True)
    receptor, reference_ligand = _stage_fixtures(output_dir)
    _write_config(output_dir / "run_config.yaml", receptor=receptor, reference_ligand=reference_ligand)
    _write_sites(output_dir / "sites.tsv")
    _write_predictions(output_dir / "predictions.tsv")
    _write_top_hits(output_dir / "top_hits.csv")
    _write_sqlite(output_dir / "results.sqlite")
    outputs = generate_report(output_dir)
    print(f"Quickstart run directory: {output_dir}")
    print(f"Report: {outputs['report_html']}")


if __name__ == "__main__":
    main()
