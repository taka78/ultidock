#!/usr/bin/env python3
"""Evaluate auto-detected UltiDock cavity centers against a crystal ligand."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


DEFAULT_HIT_THRESHOLD_A = 4.0
_SITE_ID_RE = re.compile(r"^S?(\d+)$", re.IGNORECASE)


def calculate_mol2_center(mol2_path: Path) -> tuple[float, float, float]:
    """Calculate the geometric center of a TRIPOS MOL2 ligand."""
    lines = mol2_path.read_text(encoding="utf-8").splitlines()

    in_atom_block = False
    x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
    count = 0

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_block = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom_block:
            break

        if not in_atom_block:
            continue

        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            x_sum += float(parts[2])
            y_sum += float(parts[3])
            z_sum += float(parts[4])
            count += 1
        except ValueError:
            continue

    if count == 0:
        raise ValueError(f"No atoms found in {mol2_path}")

    return (x_sum / count, y_sum / count, z_sum / count)


def _looks_like_header(parts: list[str]) -> bool:
    lowered = [part.strip().lower() for part in parts]
    return any(token in {"site", "site_id", "id", "center_x", "cx"} for token in lowered)


def _parse_center_row(parts: list[str]) -> tuple[str, tuple[float, float, float]] | None:
    if len(parts) < 5:
        return None
    site_id = parts[1].strip()
    if not site_id:
        return None
    try:
        return site_id, (float(parts[2]), float(parts[3]), float(parts[4]))
    except ValueError:
        return None


def parse_centers_tsv(tsv_path: Path) -> dict[str, tuple[float, float, float]]:
    """Parse a legacy or headered centers.tsv file into site centers."""
    sites: dict[str, tuple[float, float, float]] = {}
    header: list[str] | None = None

    for raw in tsv_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t") if "\t" in line else [part.strip() for part in line.split(",")]
        if header is None and _looks_like_header(parts):
            header = [part.strip().lower() for part in parts]
            continue

        if header is None:
            parsed = _parse_center_row(parts)
            if parsed is not None:
                site_id, center = parsed
                sites[site_id] = center
            continue

        row = {key: value for key, value in zip(header, parts)}
        site_id = row.get("site") or row.get("site_id") or row.get("id")
        if not site_id:
            continue
        try:
            cx = float(row.get("center_x") or row.get("x") or row.get("cx") or "")
            cy = float(row.get("center_y") or row.get("y") or row.get("cy") or "")
            cz = float(row.get("center_z") or row.get("z") or row.get("cz") or "")
        except ValueError:
            continue
        sites[site_id] = (cx, cy, cz)
    return sites


def euclidean_distance(p1: tuple[float, float, float], p2: tuple[float, float, float]) -> float:
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def _site_sort_key(site_id: str) -> tuple[int, str]:
    match = _SITE_ID_RE.match(site_id.strip())
    if match:
        return int(match.group(1)), site_id
    return (10**9, site_id)


def evaluate_sites(
    *,
    reference_ligand_path: Path,
    centers_tsv_path: Path,
    hit_threshold_a: float = DEFAULT_HIT_THRESHOLD_A,
) -> dict[str, Any]:
    """Return cavity-finder accuracy metrics against a crystal ligand centroid."""
    if not reference_ligand_path.exists():
        raise FileNotFoundError(f"Crystal ligand not found: {reference_ligand_path}")
    if not centers_tsv_path.exists():
        raise FileNotFoundError(f"centers.tsv not found: {centers_tsv_path}")

    true_center = calculate_mol2_center(reference_ligand_path)
    predicted_sites = parse_centers_tsv(centers_tsv_path)
    if not predicted_sites:
        raise ValueError(f"No predicted sites found in {centers_tsv_path}")

    site_rows: list[dict[str, Any]] = []
    for site_order, site_id in enumerate(sorted(predicted_sites, key=_site_sort_key), start=1):
        center = predicted_sites[site_id]
        distance = euclidean_distance(true_center, center)
        site_rows.append(
            {
                "site_id": site_id,
                "site_order": site_order,
                "distance_a": distance,
                "center": {"x": center[0], "y": center[1], "z": center[2]},
                "hit": distance <= hit_threshold_a,
            }
        )

    best_site = min(site_rows, key=lambda row: row["distance_a"])
    return {
        "reference_ligand": str(reference_ligand_path.resolve()),
        "centers_tsv": str(centers_tsv_path.resolve()),
        "reference_center": {"x": true_center[0], "y": true_center[1], "z": true_center[2]},
        "hit_threshold_a": hit_threshold_a,
        "n_sites": len(site_rows),
        "best_site_id": best_site["site_id"],
        "best_distance_a": best_site["distance_a"],
        "best_site_order": best_site["site_order"],
        "any_site_success": bool(best_site["hit"]),
        "sites": site_rows,
    }


def write_evaluation_json(path: Path, evaluation: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(evaluation, indent=2), encoding="utf-8")
    return path


def write_evaluation_csv(path: Path, evaluation: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["site_id", "site_order", "distance_a", "hit", "center_x", "center_y", "center_z"])
        for row in evaluation["sites"]:
            writer.writerow(
                [
                    row["site_id"],
                    row["site_order"],
                    f"{row['distance_a']:.4f}",
                    "1" if row["hit"] else "0",
                    f"{row['center']['x']:.4f}",
                    f"{row['center']['y']:.4f}",
                    f"{row['center']['z']:.4f}",
                ]
            )
    return path


def _default_paths(dataset_dir: Path, macro_mol_dir: Path) -> tuple[Path, Path]:
    crystal_ligand_path = dataset_dir / "crystal_ligand.mol2"
    centers_tsv_path = macro_mol_dir / "receptor" / "centers.tsv"
    return crystal_ligand_path, centers_tsv_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare auto-detected UltiDock cavity centers to a crystal ligand centroid.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dataset_dir", type=Path, help="Directory containing crystal_ligand.mol2.")
    parser.add_argument("macro_mol_dir", type=Path, help="Directory containing receptor/centers.tsv.")
    parser.add_argument(
        "--hit-threshold",
        type=float,
        default=DEFAULT_HIT_THRESHOLD_A,
        help="Distance threshold in Angstrom for a site to count as a hit.",
    )
    parser.add_argument("--output-json", type=Path, help="Optional JSON output path.")
    parser.add_argument("--output-csv", type=Path, help="Optional CSV output path.")
    args = parser.parse_args()

    crystal_ligand_path, centers_tsv_path = _default_paths(args.dataset_dir, args.macro_mol_dir)
    try:
        evaluation = evaluate_sites(
            reference_ligand_path=crystal_ligand_path,
            centers_tsv_path=centers_tsv_path,
            hit_threshold_a=args.hit_threshold,
        )
    except Exception as exc:
        print(f"[ERROR] {exc}")
        sys.exit(1)

    ref = evaluation["reference_center"]
    print(
        "True Orthosteric Center (Crystal Ligand): "
        f"({ref['x']:.2f}, {ref['y']:.2f}, {ref['z']:.2f})"
    )
    print(f"\nPredicted Sites ({evaluation['n_sites']} found):")
    for row in evaluation["sites"]:
        center = row["center"]
        marker = "★ HIT" if row["hit"] else ""
        print(
            f"  {row['site_id']}: {row['distance_a']:6.2f} A error | "
            f"Center: ({center['x']:.2f}, {center['y']:.2f}, {center['z']:.2f}) | {marker}"
        )

    print("\n" + "=" * 50)
    print("CAVITY EVALUATION SUMMARY")
    print("=" * 50)
    print(f"Distance Threshold        : <= {evaluation['hit_threshold_a']:.2f} A")
    print(
        f"Best Site                  : {evaluation['best_site_id']} "
        f"(Site order {evaluation['best_site_order']}, Error: {evaluation['best_distance_a']:.2f} A)"
    )
    print(f"Any-Site Success           : {'YES' if evaluation['any_site_success'] else 'NO'}")
    print("=" * 50)

    if args.output_json:
        write_evaluation_json(args.output_json, evaluation)
    if args.output_csv:
        write_evaluation_csv(args.output_csv, evaluation)


if __name__ == "__main__":
    main()
