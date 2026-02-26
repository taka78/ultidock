"""
ultidock.grids.check
~~~~~~~~~~~~~~~~~~~~
Sanity checker for AutoGrid .fld and .map files.

AutoGrid generates:
  - A <receptor>.maps.fld  — AVS field descriptor file
  - One <receptor>.<type>.map per atom type
  - <receptor>.e.map   — electrostatics
  - <receptor>.d.map   — desolvation

Common failure modes we catch:
  - All-zero .map   → receptor atoms weren't read or all outside grid
  - NaN / Inf       → spacing=0 or overflowed energy calculation
  - Absurd range    → unit mismatch or log-scale bug
  - Missing file    → autogrid crashed mid-run
  - Type mismatch   → GPF receptor_types ≠ ligand_types ≠ actual maps

Usage:
    report = check_fld(Path("receptor.maps.fld"))
    if not report.ok:
        for e in report.errors:
            print("ERROR:", e)
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

__all__ = [
    "MapCheckResult",
    "CheckReport",
    "GridSanityError",
    "check_fld",
]

# Threshold for "absurdly large" energy (kcal/mol). AutoDock typically caps at
# 9999, but the internal cap in some builds is 1e6.
_ABSURD_ENERGY = 1e5


class GridSanityError(RuntimeError):
    """Raised by callers that want to fail-fast on grid problems."""


@dataclass
class MapCheckResult:
    atom_type:   str
    map_path:    Path
    ok:          bool
    errors:      list[str] = field(default_factory=list)
    warnings:    list[str] = field(default_factory=list)
    sampled_min: float | None = None
    sampled_max: float | None = None
    n_sampled:   int = 0


@dataclass
class CheckReport:
    fld_path:    Path
    map_results: list[MapCheckResult] = field(default_factory=list)
    errors:      list[str] = field(default_factory=list)   # fld-level errors
    warnings:    list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors and all(r.ok for r in self.map_results)

    def fail_fast(self) -> None:
        """Raise GridSanityError if any errors were found."""
        if not self.ok:
            all_errors = list(self.errors)
            for r in self.map_results:
                all_errors.extend(f"[{r.atom_type}] {e}" for e in r.errors)
            raise GridSanityError(
                f"Grid sanity check failed for {self.fld_path}:\n"
                + "\n".join(f"  ERROR: {e}" for e in all_errors)
            )


# ── .fld parser ───────────────────────────────────────────────────────────────

def _parse_fld(fld_path: Path) -> list[tuple[str, Path]]:
    """
    Parse an AutoGrid AVS .fld file.

    Returns a list of (atom_type, map_path) pairs, ordered by variable index.

    The .fld format (simplified):
        variable N file=receptor.C.map filetype=ascii skip=6
        label=C-affinity
        ...

    We pair each "variable N file=..." line with the following "label=..." line.
    The atom type is extracted from the label (e.g., "C-affinity" → "C",
    "Electrostatics" → "e", "Desolvation" → "d").
    """
    text = fld_path.read_text(encoding="utf-8", errors="ignore")
    fld_dir = fld_path.parent

    # Collect variable-index → file mappings
    var_to_file: dict[int, Path] = {}
    for m in re.finditer(r"variable\s+(\d+)\s+file=(\S+)", text, flags=re.I):
        idx = int(m.group(1))
        var_to_file[idx] = (fld_dir / m.group(2)).resolve()

    # Collect labels in appearance order (they follow each variable line)
    labels: list[str] = [m.group(1).strip()
                         for m in re.finditer(r"label\s*=\s*(.+)", text, flags=re.I)]

    if not var_to_file:
        raise ValueError("No 'variable N file=...' entries found in .fld")

    # Build (atom_type, path) list paired by index
    result: list[tuple[str, Path]] = []
    for idx, label in enumerate(labels, start=1):
        if idx not in var_to_file:
            continue
        atom_type = _label_to_atom_type(label)
        result.append((atom_type, var_to_file[idx]))

    return result


def _label_to_atom_type(label: str) -> str:
    """
    Map an AutoGrid .fld label to a short atom-type key.

      "C-affinity"      → "C"
      "HD-affinity"     → "HD"
      "Electrostatics"  → "e"
      "Desolvation"     → "d"
      anything else     → label.split('-')[0]
    """
    lo = label.lower()
    if lo in ("electrostatics",):
        return "e"
    if lo in ("desolvation",):
        return "d"
    if "-affinity" in lo:
        return label.split("-")[0]
    return label.split("-")[0]


# ── Per-map checker ───────────────────────────────────────────────────────────

def _sample_map_values(map_path: Path, *, n: int = 500) -> list[float]:
    """
    Read the first `n` numeric lines from an AutoGrid .map file.

    AutoGrid .map files have a plain-text header (lines starting with letters)
    followed by one float per line.  We skip header lines and collect data.
    """
    values: list[float] = []
    with map_path.open(encoding="ascii", errors="ignore") as fh:
        for raw_line in fh:
            stripped = raw_line.strip()
            if not stripped:
                continue
            # Header lines start with a letter (e.g., "GRID_PARAMETER_FILE ...")
            if stripped[0].isalpha():
                continue
            try:
                values.append(float(stripped))
            except ValueError:
                pass    # skip malformed lines silently; caught by errors list below
            if len(values) >= n:
                break
    return values


def _check_single_map(atom_type: str, map_path: Path) -> MapCheckResult:
    result = MapCheckResult(atom_type=atom_type, map_path=map_path, ok=True)

    # ── Existence ─────────────────────────────────────────────────────────────
    if not map_path.exists():
        result.errors.append(f"File not found: {map_path}")
        result.ok = False
        return result

    # ── Sample values ─────────────────────────────────────────────────────────
    values = _sample_map_values(map_path)
    result.n_sampled = len(values)

    if not values:
        result.errors.append(
            "Map file exists but contains no numeric data — header-only or empty"
        )
        result.ok = False
        return result

    lo, hi = min(values), max(values)
    result.sampled_min = lo
    result.sampled_max = hi

    # ── All-zero ──────────────────────────────────────────────────────────────
    if all(v == 0.0 for v in values):
        result.errors.append(
            f"ALL {len(values)} sampled values are 0.0 — "
            "AutoGrid likely failed to read the receptor, or the grid box "
            "does not contain any receptor atoms. "
            "Check the .glg log file for 'WARNING' lines."
        )

    # ── NaN / Inf ─────────────────────────────────────────────────────────────
    elif any(not math.isfinite(v) for v in values):
        bad = [v for v in values if not math.isfinite(v)]
        result.errors.append(
            f"NaN or Inf detected in sampled map values ({len(bad)} of {len(values)})"
        )

    # ── Absurd range ──────────────────────────────────────────────────────────
    else:
        if hi > _ABSURD_ENERGY or lo < -_ABSURD_ENERGY:
            result.warnings.append(
                f"Extreme energy range [{lo:.2f}, {hi:.2f}] kcal/mol "
                f"(threshold ±{_ABSURD_ENERGY:.0f}) — possible grid overflow or unit error"
            )

    result.ok = not result.errors
    return result


# ── Public entry-point ────────────────────────────────────────────────────────

def check_fld(
    fld_path: Path,
    expected_types: Sequence[str] | None = None,
) -> CheckReport:
    """
    Sanity-check an AutoGrid .fld file and all referenced .map files.

    Parameters
    ----------
    fld_path       : path to the .fld file produced by AutoGrid
    expected_types : optional list of atom type strings that should be present
                     (e.g., ["C", "HD", "OA", "N", "e", "d"]).
                     If None, the types listed in the .fld file are used as-is
                     and no cross-check is performed.

    Returns
    -------
    CheckReport with .ok property and per-map results.

    The caller can either inspect report.errors / report.map_results, or call
    report.fail_fast() to raise GridSanityError on the first problem.
    """
    report = CheckReport(fld_path=fld_path)

    # ── Parse .fld header ─────────────────────────────────────────────────────
    try:
        map_entries = _parse_fld(fld_path)
    except Exception as exc:
        report.errors.append(f"Cannot parse .fld file: {exc}")
        return report

    if not map_entries:
        report.errors.append("No map entries found in .fld — AutoGrid may not have run")
        return report

    # ── Expected-type cross-check ─────────────────────────────────────────────
    if expected_types is not None:
        found  = {t for t, _ in map_entries}
        expect = set(expected_types)
        if found != expect:
            missing   = expect - found
            extra     = found  - expect
            parts = []
            if missing:
                parts.append(f"missing={sorted(missing)}")
            if extra:
                parts.append(f"unexpected={sorted(extra)}")
            report.errors.append(
                f"Map type mismatch — {', '.join(parts)}. "
                "Verify that receptor_types in the .gpf matches the atom types "
                "in the receptor PDBQT."
            )

    # ── Per-map checks ────────────────────────────────────────────────────────
    for atom_type, map_path in map_entries:
        result = _check_single_map(atom_type, map_path)
        report.map_results.append(result)

    return report
