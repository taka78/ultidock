#!/usr/bin/env python3
"""Parse AutoDock-GPU log files and summarize timing and memory metrics."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict, fields
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from db_manager import DockingDatabaseManager


CUDA_SETUP_RE = re.compile(r"CUDA Setup time\s+([0-9]*\.?[0-9]+)s")
REST_SETUP_RE = re.compile(r"Rest of Setup time\s+([0-9]*\.?[0-9]+)s")
SHUTDOWN_RE = re.compile(r"Shutdown time\s+([0-9]*\.?[0-9]+)s")
JOB_TIMING_RE = re.compile(
    r"Job #(?P<job>\d+) took\s+([0-9]*\.?[0-9]+)\s+sec after waiting\s+([0-9]*\.?[0-9]+)\s+sec for setup"
)
RUN_TIME_RE = re.compile(r"Run time of entire job set.*:\s+([0-9]*\.?[0-9]+)\s+sec")
PROCESSING_RE = re.compile(r"Processing time:\s+([0-9]*\.?[0-9]+)\s+sec")
MEMORY_RE = re.compile(r"Memory usage:\s+([0-9]*\.?[0-9]+)\s+MB")
USING_GRID_RE = re.compile(r"Using grid:\s*(.*)")
PROCESSING_LIGAND_RE = re.compile(r"Processing\s+(.*\.pdbqt)")
STDERR_LIGAND_RE = re.compile(r"(?P<ligand>(?:[A-Za-z]:)?[^\s]+\.pdbqt)$")


METRIC_COLUMN_MAP = {
    "cuda_setup_time_s": "cuda_setup_s",
    "rest_setup_time_s": "rest_setup_s",
    "total_setup_time_s": "total_setup_s",
    "job_wait_time_s": "job_wait_s",
    "job_duration_s": "job_duration_s",
    "run_time_s": "run_time_s",
    "processing_time_s": "processing_time_s",
    "docking_plus_processing_s": "run_plus_processing_s",
    "shutdown_time_s": "shutdown_s",
    "memory_usage_mb": "gpu_memory_mb",
}

METRIC_DB_COLUMNS = list(METRIC_COLUMN_MAP.values())


def normalize_ligand_name(path_str: Optional[str]) -> Optional[str]:
    """Return the ligand identifier prefix (sans extension) for DB matching."""
    if not path_str:
        return None

    try:
        name = Path(path_str).name
    except Exception:
        name = str(path_str)

    if name.lower().endswith(".pdbqt"):
        name = name[:-6]

    return name or None


def ligand_prefix_from_db_name(name: Optional[str]) -> str:
    if not name:
        return ""

    if "-Model" in name:
        return name.split("-Model", 1)[0]

    if name.lower().endswith(".pdbqt"):
        return name[:-6]

    return name


@dataclass
class DockingRecord:
    grid_path: Optional[str] = None
    ligand_path: Optional[str] = None
    ligand_name: Optional[str] = None
    cuda_setup_time_s: Optional[float] = None
    rest_setup_time_s: Optional[float] = None
    shutdown_time_s: Optional[float] = None
    job_duration_s: Optional[float] = None
    job_wait_time_s: Optional[float] = None
    run_time_s: Optional[float] = None
    processing_time_s: Optional[float] = None
    memory_usage_mb: Optional[float] = None
    job_index: Optional[int] = None
    site: Optional[str] = field(default=None, init=False)
    binding_site: Optional[str] = field(default=None, init=False)

    def finalize(self) -> None:
        if self.grid_path:
            path = Path(self.grid_path)
            # prefer the parent folder (e.g. S1, S2) if it exists
            parent_name = path.parent.name
            self.site = parent_name or path.name
        else:
            self.site = None

        if self.site:
            match = re.search(r"[Ss]?(\d+)", self.site)
            if match:
                self.binding_site = match.group(1)
        elif self.grid_path:
            match = re.search(r"/(?:S|s)(\d+)/", self.grid_path)
            if match:
                self.binding_site = match.group(1)

    @property
    def total_setup_time_s(self) -> Optional[float]:
        if self.cuda_setup_time_s is None and self.rest_setup_time_s is None:
            return None
        total = 0.0
        if self.cuda_setup_time_s is not None:
            total += self.cuda_setup_time_s
        if self.rest_setup_time_s is not None:
            total += self.rest_setup_time_s
        return total

    @property
    def docking_plus_processing_s(self) -> Optional[float]:
        if self.run_time_s is None and self.processing_time_s is None:
            return None
        total = 0.0
        if self.run_time_s is not None:
            total += self.run_time_s
        if self.processing_time_s is not None:
            total += self.processing_time_s
        return total

    def as_serializable_dict(self) -> dict:
        result = asdict(self)
        result["total_setup_time_s"] = self.total_setup_time_s
        result["docking_plus_processing_s"] = self.docking_plus_processing_s
        result["site"] = self.site
        result["binding_site"] = self.binding_site
        return result

    def as_db_update(self) -> Dict[str, float]:
        payload = {}
        for attr, column in METRIC_COLUMN_MAP.items():
            value = getattr(self, attr, None)
            if value is not None:
                payload[column] = value
        return payload


def parse_records(lines: Iterable[str]) -> List[DockingRecord]:
    records: List[DockingRecord] = []
    current: Optional[DockingRecord] = None
    current_section: Optional[str] = None
    last_grid: Optional[str] = None
    last_ligand_path: Optional[str] = None
    last_ligand_name: Optional[str] = None

    for raw_line in lines:
        line = raw_line.strip()

        grid_match = USING_GRID_RE.match(line)
        if grid_match:
            last_grid = grid_match.group(1).strip() or None

        ligand_match = PROCESSING_LIGAND_RE.match(line)
        if ligand_match:
            last_ligand_path = ligand_match.group(1).strip() or None
            last_ligand_name = normalize_ligand_name(last_ligand_path)

        stderr_ligand_match = STDERR_LIGAND_RE.match(line)
        if stderr_ligand_match:
            last_ligand_path = stderr_ligand_match.group("ligand").strip() or None
            last_ligand_name = normalize_ligand_name(last_ligand_path)

        if line == "[AutoDock-GPU stdout]":
            current_section = "stdout"
            current = DockingRecord(
                grid_path=last_grid,
                ligand_path=last_ligand_path,
                ligand_name=last_ligand_name,
            )
            records.append(current)
            continue
        if line == "[AutoDock-GPU stderr]":
            current_section = "stderr"
            continue
        if line.startswith("[") and line.endswith("]") and "AutoDock-GPU" not in line:
            current_section = None
            current = None
            continue

        if current is None or current_section is None:
            continue

        if current_section == "stdout":
            match = CUDA_SETUP_RE.search(line)
            if match:
                current.cuda_setup_time_s = float(match.group(1))
                continue

            match = REST_SETUP_RE.search(line)
            if match:
                current.rest_setup_time_s = float(match.group(1))
                continue

            match = SHUTDOWN_RE.search(line)
            if match:
                current.shutdown_time_s = float(match.group(1))
                continue

            match = JOB_TIMING_RE.search(line)
            if match:
                current.job_index = int(match.group("job"))
                current.job_duration_s = float(match.group(2))
                current.job_wait_time_s = float(match.group(3))
                continue

            match = RUN_TIME_RE.search(line)
            if match:
                current.run_time_s = float(match.group(1))
                continue

            match = PROCESSING_RE.search(line)
            if match:
                current.processing_time_s = float(match.group(1))
                continue

        elif current_section == "stderr":
            match = MEMORY_RE.search(line)
            if match:
                current.memory_usage_mb = float(match.group(1))
                continue

    for record in records:
        record.finalize()

    return records


def format_float(value: Optional[float]) -> str:
    return f"{value:.3f}" if value is not None else "-"


def short_name(path_str: Optional[str]) -> str:
    if not path_str:
        return "-"
    try:
        return Path(path_str).name
    except Exception:
        return str(path_str)


def format_table(records: List[DockingRecord]) -> str:
    headers = [
        "#",
        "ligand_name",
        "site",
        "binding",
        "cuda_setup_s",
        "rest_setup_s",
        "total_setup_s",
        "wait_s",
        "job_s",
        "run_s",
        "processing_s",
        "run+proc_s",
        "shutdown_s",
        "mem_mb",
    ]

    rows = []
    for idx, record in enumerate(records, start=1):
        display = record.ligand_name or short_name(record.ligand_path)
        rows.append([
            str(idx),
            display,
            record.site or "-",
            record.binding_site or "-",
            format_float(record.cuda_setup_time_s),
            format_float(record.rest_setup_time_s),
            format_float(record.total_setup_time_s),
            format_float(record.job_wait_time_s),
            format_float(record.job_duration_s),
            format_float(record.run_time_s),
            format_float(record.processing_time_s),
            format_float(record.docking_plus_processing_s),
            format_float(record.shutdown_time_s),
            format_float(record.memory_usage_mb),
        ])

    widths = [len(header) for header in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def join_row(row: List[str]) -> str:
        return "  ".join(cell.ljust(widths[idx]) for idx, cell in enumerate(row))

    if not rows:
        return "No AutoDock-GPU runs found in log."

    lines = [join_row(headers)]
    lines.append(
        join_row(["-" * widths[idx] for idx, _ in enumerate(headers)])
    )
    for row in rows:
        lines.append(join_row(row))

    return "\n".join(lines)


def summarize(records: List[DockingRecord]) -> List[str]:
    summary_lines: List[str] = []
    if not records:
        summary_lines.append("No AutoDock-GPU runs found in log.")
        return summary_lines

    total_run = sum((r.run_time_s or 0.0) for r in records)
    total_processing = sum((r.processing_time_s or 0.0) for r in records)
    total_run_proc = sum((r.docking_plus_processing_s or 0.0) for r in records)
    total_setup = sum((r.total_setup_time_s or 0.0) for r in records)

    mem_values = [r.memory_usage_mb for r in records if r.memory_usage_mb is not None]

    summary_lines.append(f"Entries parsed: {len(records)}")
    unique_keys = {
        (r.ligand_name or "", r.binding_site or "")
        for r in records
        if r.ligand_name
    }
    if unique_keys:
        summary_lines.append(f"Unique ligand/site combinations: {len(unique_keys)}")
    summary_lines.append(f"Total CUDA+rest setup time: {total_setup:.3f} s")
    summary_lines.append(f"Total run time: {total_run:.3f} s")
    summary_lines.append(f"Total processing time: {total_processing:.3f} s")
    summary_lines.append(f"Total docking time incl. processing: {total_run_proc:.3f} s")

    if mem_values:
        avg_mem = sum(mem_values) / len(mem_values)
        max_mem = max(mem_values)
        summary_lines.append(f"Average memory usage: {avg_mem:.2f} MB")
        summary_lines.append(f"Peak memory usage: {max_mem:.2f} MB")

    return summary_lines


def merge_records(target: DockingRecord, incoming: DockingRecord) -> DockingRecord:
    if not target.ligand_path and incoming.ligand_path:
        target.ligand_path = incoming.ligand_path
    if not target.ligand_name and incoming.ligand_name:
        target.ligand_name = incoming.ligand_name
    if not target.grid_path and incoming.grid_path:
        target.grid_path = incoming.grid_path
    if incoming.site and not target.site:
        target.site = incoming.site
    if incoming.binding_site and not target.binding_site:
        target.binding_site = incoming.binding_site

    writable_attrs = _DOCKING_RECORD_WRITABLE_ATTRS
    for attr in METRIC_COLUMN_MAP.keys():
        if attr not in writable_attrs:
            continue
        value = getattr(incoming, attr, None)
        if value is not None:
            setattr(target, attr, value)

    if incoming.job_index is not None:
        target.job_index = incoming.job_index

    return target


def aggregate_by_ligand(records: List[DockingRecord]) -> Tuple[Dict[Tuple[str, Optional[str]], DockingRecord], Dict[Tuple[str, Optional[str]], int]]:
    aggregated: Dict[Tuple[str, Optional[str]], DockingRecord] = {}
    duplicates: Dict[Tuple[str, Optional[str]], int] = defaultdict(int)

    for record in records:
        if not record.ligand_name:
            continue
        key = (record.ligand_name, record.binding_site)
        if key in aggregated:
            duplicates[key] += 1
            aggregated[key] = merge_records(aggregated[key], record)
        else:
            aggregated[key] = record

    return aggregated, duplicates


_DOCKING_RECORD_WRITABLE_ATTRS = {field.name for field in fields(DockingRecord)}


def update_database_from_records(
    records: List[DockingRecord],
    db_filename: Optional[str],
    dry_run: bool = False,
) -> Dict[str, object]:
    aggregated, duplicates = aggregate_by_ligand(records)

    manager = DockingDatabaseManager(db_filename or 'ultidock_results.db')
    stats = {
        "matched_keys": 0,
        "updated_rows": 0,
        "missing": [],
        "skipped": [],
        "duplicate_keys": duplicates,
        "touched_rows": [],
    }

    try:
        for (ligand_prefix, binding_site), record in aggregated.items():
            metrics = record.as_db_update()
            if not metrics:
                stats["skipped"].append((ligand_prefix, binding_site, "no metrics"))
                continue

            matches = manager.fetch_rows_by_ligand_prefix(ligand_prefix, binding_site)
            match_count = len(matches)

            if match_count == 0:
                stats["missing"].append((ligand_prefix, binding_site))
                continue

            stats["matched_keys"] += 1
            if dry_run:
                stats["updated_rows"] += match_count
                stats["touched_rows"].extend(matches)
                continue

            updated = manager.update_metrics_by_ligand_prefix(
                ligand_prefix, binding_site, metrics
            )
            stats["updated_rows"] += updated

            refreshed = manager.fetch_rows_by_ligand_prefix(ligand_prefix, binding_site)
            stats["touched_rows"].extend(refreshed or matches)
    finally:
        manager.close()

    return stats


def _aggregate_rows_for_csv(
    rows: List[Dict[str, object]],
    key_fn,
    label: str,
) -> List[Dict[str, object]]:
    grouped: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[key_fn(row)].append(row)

    summaries: List[Dict[str, object]] = []
    for key, group_rows in grouped.items():
        entry: Dict[str, object] = {label: key, "entry_count": len(group_rows)}
        for column in METRIC_DB_COLUMNS:
            values = [r.get(column) for r in group_rows if r.get(column) is not None]
            if values:
                total = sum(values)
                entry[f"{column}_avg"] = total / len(values)
                entry[f"{column}_sum"] = total
            else:
                entry[f"{column}_avg"] = None
                entry[f"{column}_sum"] = None
        summaries.append(entry)

    summaries.sort(key=lambda item: item.get(label) or "")
    return summaries


def write_csv_summaries(rows: List[Dict[str, object]], log_path: Path) -> List[Path]:
    if not rows:
        return []

    output_dir = log_path.parent
    base = log_path.stem

    generated: List[Path] = []

    detail_fields = ["id", "ligand_name", "binding_site", *METRIC_DB_COLUMNS]
    detail_path = output_dir / f"{base}_docking_metrics.csv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields)
        writer.writeheader()
        for row in sorted(
            rows,
            key=lambda r: ((r.get("binding_site") or ""), r.get("ligand_name") or ""),
        ):
            writer.writerow({field: row.get(field) for field in detail_fields})
    generated.append(detail_path)

    site_rows = _aggregate_rows_for_csv(
        rows,
        key_fn=lambda r: r.get("binding_site") or "",
        label="binding_site",
    )
    site_fields = [
        "binding_site",
        "entry_count",
        *[f"{column}_avg" for column in METRIC_DB_COLUMNS],
        *[f"{column}_sum" for column in METRIC_DB_COLUMNS],
    ]
    site_path = output_dir / f"{base}_site_summary.csv"
    with site_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=site_fields)
        writer.writeheader()
        for row in site_rows:
            writer.writerow({field: row.get(field) for field in site_fields})
    generated.append(site_path)

    ligand_rows = _aggregate_rows_for_csv(
        rows,
        key_fn=lambda r: ligand_prefix_from_db_name(r.get("ligand_name")),
        label="ligand_prefix",
    )
    ligand_fields = [
        "ligand_prefix",
        "entry_count",
        *[f"{column}_avg" for column in METRIC_DB_COLUMNS],
        *[f"{column}_sum" for column in METRIC_DB_COLUMNS],
    ]
    ligand_path = output_dir / f"{base}_ligand_summary.csv"
    with ligand_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ligand_fields)
        writer.writeheader()
        for row in ligand_rows:
            writer.writerow({field: row.get(field) for field in ligand_fields})
    generated.append(ligand_path)

    return generated


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse AutoDock-GPU logs and summarize timing / memory metrics."
    )
    parser.add_argument("logfile", help="Path to the log file to parse")
    parser.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="Emit results as JSON instead of a formatted table",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output (implies --json)",
    )
    parser.add_argument(
        "--db",
        help="SQLite database filename or absolute path to update with parsed metrics",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and show matches without writing to the database",
    )
    parser.add_argument(
        "--no-table",
        action="store_true",
        help="Skip the formatted table output (summary and DB update still run)",
    )

    args = parser.parse_args()

    path = Path(args.logfile)
    if not path.is_file():
        parser.error(f"Log file not found: {path}")

    records = parse_records(path.read_text().splitlines())

    if args.pretty:
        args.as_json = True

    if args.as_json:
        payload = [record.as_serializable_dict() for record in records]
        json_kwargs = {"indent": 2, "sort_keys": True} if args.pretty else {}
        json.dump(payload, sys.stdout, **json_kwargs)
        print()
        return

    if not records:
        print("No AutoDock-GPU runs found in log.")
        if args.db:
            print("No database updates were attempted because no runs were detected.")
        return

    if not args.no_table:
        table = format_table(records)
        print(table)
        print()

    for line in summarize(records):
        print(line)

    if args.db:
        print()
        stats = update_database_from_records(records, args.db, dry_run=args.dry_run)

        duplicates_total = sum(stats["duplicate_keys"].values()) if stats["duplicate_keys"] else 0
        if args.dry_run:
            print(
                f"[DRY-RUN] {stats['matched_keys']} ligand/site combinations would be updated in"
                f" {stats['updated_rows']} rows."
            )
        else:
            print(
                f"Database update complete: {stats['updated_rows']} rows touched across"
                f" {stats['matched_keys']} ligand/site combinations."
            )

        if duplicates_total:
            print(
                f"Note: detected {duplicates_total} additional log entries for"
                " existing ligand/site combinations (latest values were used)."
            )

        if stats["missing"]:
            preview = ", ".join(
                f"{short_name(l)} (site {s or '-'})" for l, s in stats["missing"][:5]
            )
            more = "" if len(stats["missing"]) <= 5 else " …"
            print(
                f"Warning: {len(stats['missing'])} combinations were not found in the database:"
                f" {preview}{more}"
            )

        if stats["skipped"]:
            print(
                f"Skipped {len(stats['skipped'])} combinations without usable metrics entries."
            )

        if not args.dry_run and stats["touched_rows"]:
            csv_paths = write_csv_summaries(stats["touched_rows"], path)
            if csv_paths:
                print()
                for csv_path in csv_paths:
                    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()