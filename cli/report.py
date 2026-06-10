"""Report and visualization helpers for Ultidock run directories."""

from __future__ import annotations

import csv
import html
import json
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _repo_root() -> Path:
    for candidate in [Path(__file__).resolve().parent.parent, Path.cwd(), *Path.cwd().parents]:
        if (candidate / "docking").is_dir() and (candidate / "README.md").is_file():
            return candidate
    return Path.cwd()


def _git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_repo_root(),
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def _read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _read_key_values(path: Path | None) -> dict[str, str]:
    if path is None or not path.is_file():
        return {}
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip().strip("'\"")
    return values


def _resolve_path(raw_path: str, *, base: Path) -> Path | None:
    if not raw_path:
        return None
    path = Path(raw_path)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _pdb_summary(path: Path | None) -> dict[str, Any]:
    if path is None or not path.is_file():
        return {}
    atoms = 0
    chains: set[str] = set()
    residues: set[tuple[str, str, str]] = set()
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        atoms += 1
        chain = raw[21].strip() or "_"
        resname = raw[17:20].strip()
        resseq = raw[22:26].strip()
        chains.add(chain)
        residues.add((chain, resseq, resname))
    return {
        "path": str(path),
        "atom_count": atoms,
        "chain_count": len(chains),
        "chains": ",".join(sorted(chains)),
        "residue_count": len(residues),
    }


def _csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.4f}"
        return ""
    return str(value)


def _parse_centers_tsv(path: Path) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split("\t")
        if line.startswith("#"):
            normalized = [part.strip().lstrip("#").strip() for part in parts]
            if {"site_id", "cx", "cy", "cz"}.issubset(set(normalized)):
                header = normalized
            continue
        if header is None:
            if len(parts) >= 11:
                header = [
                    "receptor",
                    "site_id",
                    "cx",
                    "cy",
                    "cz",
                    "nx",
                    "ny",
                    "nz",
                    "spacing",
                    "r_peak",
                    "F",
                    "raw_F",
                    "family",
                    "portfolio_role",
                    "selection_score",
                    "center_closeness",
                ][: len(parts)]
            else:
                continue
        row = {key: value for key, value in zip(header, parts)}
        rows.append(row)
    return rows


def _find_first(run_dir: Path, names: tuple[str, ...]) -> Path | None:
    for name in names:
        path = run_dir / name
        if path.is_file():
            return path
    for name in names:
        matches = sorted(run_dir.rglob(name))
        if matches:
            return matches[0]
    return None


def _site_rows(run_dir: Path) -> list[dict[str, str]]:
    centers_path = _find_first(run_dir, ("sites.tsv", "centers.tsv"))
    if centers_path:
        rows = _parse_centers_tsv(centers_path)
        if rows:
            return rows
    predictions_path = _find_first(run_dir, ("predictions.tsv",))
    if predictions_path:
        rows = _read_table(predictions_path)
        return [
            {
                "site_id": row.get("site_id", ""),
                "cx": row.get("center_x", ""),
                "cy": row.get("center_y", ""),
                "cz": row.get("center_z", ""),
                "F": row.get("score", ""),
                "family": row.get("method", ""),
                "portfolio_role": "prediction",
            }
            for row in rows
        ]
    return []


def _top_hits(run_dir: Path) -> list[dict[str, str]]:
    path = _find_first(run_dir, ("top_hits.csv", "scores.csv"))
    if not path:
        return []
    try:
        return _read_table(path)[:20]
    except Exception:
        return []


def _write_centers_pdb(path: Path, sites: list[dict[str, str]]) -> None:
    lines = ["REMARK CaV-EMPS predicted cavity/site centers"]
    for index, site in enumerate(sites, start=1):
        try:
            x = float(site.get("cx") or site.get("center_x") or 0.0)
            y = float(site.get("cy") or site.get("center_y") or 0.0)
            z = float(site.get("cz") or site.get("center_z") or 0.0)
        except ValueError:
            continue
        site_id = (site.get("site_id") or f"S{index}")[:4]
        lines.append(
            f"HETATM{index:5d}  C   CVE A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  "
            f"  REMARK {site_id}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_pymol(path: Path, sites: list[dict[str, str]], centers_name: str) -> None:
    lines = [
        f"load {centers_name}, cavemps_centers",
        "hide everything, cavemps_centers",
        "show spheres, cavemps_centers",
        "set sphere_scale, 1.2, cavemps_centers",
        "color cyan, cavemps_centers",
    ]
    for index, site in enumerate(sites, start=1):
        site_id = site.get("site_id") or f"S{index}"
        lines.append(f"# {site_id}: center=({site.get('cx')}, {site.get('cy')}, {site.get('cz')})")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_chimerax(path: Path, centers_name: str) -> None:
    lines = [
        f"open {centers_name}",
        "style sphere",
        "color cyan",
        "size stickRadius 0.4",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_bild(path: Path, sites: list[dict[str, str]]) -> None:
    lines = [".comment CaV-EMPS predicted cavity/site centers", ".color cyan"]
    for site in sites:
        try:
            x = float(site.get("cx") or site.get("center_x") or 0.0)
            y = float(site.get("cy") or site.get("center_y") or 0.0)
            z = float(site.get("cz") or site.get("center_z") or 0.0)
        except ValueError:
            continue
        lines.append(f".sphere {x:.3f} {y:.3f} {z:.3f} 1.2")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_visual_exports(run_dir: Path, sites: list[dict[str, str]]) -> dict[str, Path]:
    exports = {
        "cavity_centers_pdb": run_dir / "cavity_centers.pdb",
        "cavemps_sites_pml": run_dir / "cavemps_sites.pml",
        "site_boxes_pml": run_dir / "site_boxes.pml",
        "top_poses_pml": run_dir / "top_poses.pml",
        "chimerax_cxc": run_dir / "cavemps_sites.cxc",
        "cavity_centers_bild": run_dir / "cavity_centers.bild",
    }
    _write_centers_pdb(exports["cavity_centers_pdb"], sites)
    _write_pymol(exports["cavemps_sites_pml"], sites, "cavity_centers.pdb")
    _write_pymol(exports["site_boxes_pml"], sites, "cavity_centers.pdb")
    exports["top_poses_pml"].write_text(
        "# Load top docking poses here when pose files are present.\n",
        encoding="utf-8",
    )
    _write_chimerax(exports["chimerax_cxc"], "cavity_centers.pdb")
    _write_bild(exports["cavity_centers_bild"], sites)
    return exports


def _markdown_table(rows: list[dict[str, Any]], columns: list[str], *, limit: int = 20) -> str:
    if not rows:
        return "_None found._\n"
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(_csv_value(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines) + "\n"


def _markdown_to_html(markdown: str) -> str:
    body = []
    in_pre = False
    for line in markdown.splitlines():
        if line.startswith("```"):
            body.append("</pre>" if in_pre else "<pre>")
            in_pre = not in_pre
        elif in_pre:
            body.append(html.escape(line))
        elif line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("| "):
            body.append(f"<p><code>{html.escape(line)}</code></p>")
        elif line.strip():
            body.append(f"<p>{html.escape(line)}</p>")
        else:
            body.append("")
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Ultidock Report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:960px;margin:2rem auto;"
        "line-height:1.45}code,pre{background:#f6f6f6;padding:.15rem .25rem}"
        "pre{padding:1rem;overflow:auto}</style></head><body>"
        + "\n".join(body)
        + "</body></html>\n"
    )


def generate_report(run_dir: Path) -> dict[str, Path]:
    """Generate Markdown/HTML reports and visualization helper files."""

    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    sites = _site_rows(run_dir)
    top_hits = _top_hits(run_dir)
    exports = write_visual_exports(run_dir, sites)
    config_path = _find_first(run_dir, ("run_config.yaml", "run_metadata.json", "summary.json"))
    config = _read_key_values(config_path)
    receptor_summary = _pdb_summary(_resolve_path(config.get("receptor_input", ""), base=run_dir))

    markdown = [
        "# Ultidock Run Report",
        "",
        "## Summary",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Run directory: `{run_dir}`",
        f"- Git commit: `{_git_commit()}`",
        "- Site method: `cav-emps` (CaV-EMPS: Cavity detection via Electrostatic Map Pocket Scoring)",
        "- CaV-EMPS scores are ranking scores, not binding affinities.",
        "",
        "## Reproducibility",
        "",
        f"- Config/metadata: `{config_path.name if config_path else 'not found'}`",
        f"- Command: `{os.environ.get('ULTIDOCK_COMMAND', 'not recorded')}`",
        "",
        "## Receptor",
        "",
        _markdown_table([receptor_summary] if receptor_summary else [], ["path", "atom_count", "chain_count", "chains", "residue_count"], limit=1),
        "",
        "## CaV-EMPS Sites",
        "",
        _markdown_table(
            sites,
            ["site_id", "cx", "cy", "cz", "F", "raw_F", "family", "portfolio_role"],
            limit=30,
        ),
        "",
        "## Top Ligands",
        "",
        _markdown_table(top_hits, list(top_hits[0].keys())[:8] if top_hits else [], limit=20),
        "",
        "## Visualization Files",
        "",
    ]
    for label, path in exports.items():
        markdown.append(f"- `{path.name}`")
    markdown.append("")

    report_md = run_dir / "report.md"
    report_html = run_dir / "report.html"
    report_md.write_text("\n".join(markdown), encoding="utf-8")
    report_html.write_text(_markdown_to_html(report_md.read_text(encoding="utf-8")), encoding="utf-8")
    exports["report_md"] = report_md
    exports["report_html"] = report_html
    return exports
