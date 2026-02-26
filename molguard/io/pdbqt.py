"""
ultidock.io.pdbqt
~~~~~~~~~~~~~~~~~
PDBQT file linting, in-place numeric normalisation (ligand), and
deterministic receptor canonicalization.

v1 contract for LIGANDS
-----------------------
  - pdbqt_check()     — structural lint only; returns LintReport (no writes)
  - pdbqt_normalize() — reformats ATOM/HETATM numeric columns in-place
                        NEVER touches ROOT / BRANCH / ENDBRANCH / TORSDOF lines
                        NEVER reorders or renumbers atoms

v1 contract for RECEPTORS
--------------------------
  - canonicalize_receptor() — sort atoms by chemistry key, renumber serials,
                              reformat all numeric columns, strip non-essential REMARKs

PDBQT column map (0-indexed Python slices)
------------------------------------------
  line[0:6]   record type   ("ATOM  " / "HETATM")
  line[6:11]  serial        int, right-justified
  line[12:16] atom name     4 chars, left-padded
  line[16]    alt_loc       single char (' ' if none)
  line[17:20] res_name      3 chars
  line[21]    chain_id      1 char
  line[22:26] res_seq       int, right-justified
  line[30:38] X coord       8 chars, 3 decimals  ← fixedfmt zone
  line[38:46] Y coord       8 chars, 3 decimals  ← fixedfmt zone
  line[46:54] Z coord       8 chars, 3 decimals  ← fixedfmt zone
  line[54:60] occupancy     6 chars, 2 decimals  ← fixedfmt zone
  line[60:66] b_factor      6 chars, 2 decimals  ← fixedfmt zone
  line[70:76] charge        6 chars, 3 decimals  ← fixedfmt zone
  line[77:79] AD atom type  2 chars, right-justified
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from molguard.io.fixedfmt import (
    FixedFmtError,
    fmt_bfactor,
    fmt_charge,
    fmt_coord,
    fmt_occupancy,
)

__all__ = [
    "LintIssue",
    "LintReport",
    "LintError",
    "pdbqt_check",
    "pdbqt_normalize",
    "ReceptorAtom",
    "canonicalize_receptor",
]

# ── Record types that carry numeric fields we must check / reformat ───────────
_ATOM_RECORDS = frozenset({"ATOM", "HETATM"})

# ── Torsion-tree records that ligand normaliser must pass through verbatim ────
_TORSION_RECORDS = frozenset({"ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"})


# ─────────────────────────────────────────────────────────────────────────────
# Lint data structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class LintIssue:
    line_no: int     # 1-indexed
    column:  str     # e.g. "X", "charge", "atom_type", "file"
    raw:     str     # the raw field text that triggered the issue
    code:    str     # e.g. "EXPONENT", "NO_DECIMAL", "NAN_INF", "NON_ASCII"
    message: str


@dataclass
class LintReport:
    path:     Path
    errors:   list[LintIssue] = field(default_factory=list)
    warnings: list[LintIssue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


class LintError(RuntimeError):
    """Raised by pdbqt_normalize / canonicalize_receptor when lint fails."""


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _check_numeric_field(
    report: LintReport,
    line_no: int,
    line: str,
    col_name: str,
    sl: slice,
) -> None:
    """Check a single fixed-width numeric field in a PDBQT line."""
    expected_width = sl.stop - sl.start

    # Pad or truncate defensively so slices always work
    raw = line[sl.start : sl.stop] if len(line) >= sl.stop else line[sl.start:]

    if len(raw) != expected_width:
        report.errors.append(LintIssue(
            line_no, col_name, raw, "WIDTH",
            f"field width {len(raw)}, expected {expected_width}",
        ))
        return  # remaining checks do not apply to a truncated field

    stripped = raw.strip()

    # 1. Exponent notation — the top-level parser in AutoGrid sees 'e' and either
    #    reads garbage or silently truncates the mantissa.
    if "e" in stripped.lower():
        report.errors.append(LintIssue(
            line_no, col_name, raw, "EXPONENT",
            "scientific notation in numeric field",
        ))

    # 2. Missing decimal point — Fortran integer parse path kicks in silently.
    if stripped and "." not in stripped:
        report.errors.append(LintIssue(
            line_no, col_name, raw, "NO_DECIMAL",
            "no decimal point in numeric field",
        ))

    # 3. NaN / Inf literals (text-level before float conversion)
    if stripped.lower() in {"nan", "inf", "-inf", "+inf", "infinity", "-infinity"}:
        report.errors.append(LintIssue(
            line_no, col_name, raw, "NAN_INF",
            f"non-finite literal: {stripped!r}",
        ))
        return

    # 4. Parse check — catches "****" overflow markers and other unparseable forms
    if stripped:
        try:
            v = float(stripped)
        except ValueError:
            report.errors.append(LintIssue(
                line_no, col_name, raw, "PARSE_FAIL",
                f"cannot parse as float: {stripped!r}",
            ))
            return

        # 5. isfinite guard (catches edge cases like float("nan") from numpy writes)
        if not math.isfinite(v):
            report.errors.append(LintIssue(
                line_no, col_name, raw, "NAN_INF",
                f"parsed value is non-finite: {v}",
            ))


# Column slices for ATOM/HETATM fields (0-indexed, half-open: line[start:stop])
_COORD_FIELDS: list[tuple[str, slice]] = [
    ("X",          slice(30, 38)),
    ("Y",          slice(38, 46)),
    ("Z",          slice(46, 54)),
    ("occupancy",  slice(54, 60)),
    ("b_factor",   slice(60, 66)),
    ("charge",     slice(70, 76)),
]

# B-factor threshold above which the field fills completely (no leading space),
# causing whitespace-tokenizing parsers to merge occ+bfac into one token.
_BFACTOR_WARN_LIMIT = 99.99



# ─────────────────────────────────────────────────────────────────────────────
# Part B — Ligand check & normalize (torsion tree untouched)
# ─────────────────────────────────────────────────────────────────────────────

def pdbqt_check(path: Path) -> LintReport:
    """
    Lint a PDBQT file (ligand or receptor) and return a structured LintReport.

    Checks performed on every ATOM/HETATM line:
      - ASCII-only file (file-level, checked first)
      - Exponent notation in numeric columns
      - Missing decimal point in numeric columns
      - NaN/Inf literals and parsed non-finite values
      - Unparseable float fields ("****", etc.)
      - Missing AD atom type (cols 77-79)
      - Short lines (warning, not error)

    This function never writes anything.  Errors are collected across all lines
    before returning (fail-slow lint).
    """
    report = LintReport(path=path)

    # ── ASCII gate ─────────────────────────────────────────────────────────────
    try:
        raw_bytes = path.read_bytes()
        text_lines = raw_bytes.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        report.errors.append(LintIssue(0, "file", "", "NON_ASCII", str(exc)))
        return report   # cannot continue — byte positions are unknown

    for i, line in enumerate(text_lines, start=1):
        rec = line[:6].strip()
        if rec not in _ATOM_RECORDS:
            continue

        # ── Numeric fields ─────────────────────────────────────────────────────
        for col_name, sl in _COORD_FIELDS:
            _check_numeric_field(report, i, line, col_name, sl)

        # ── B-factor overflow warning ──────────────────────────────────────────
        # B-factor > 99.99 fills all 6 chars of the field with no leading space.
        # Autogrid4's whitespace tokenizer then merges occ+bfac into one token,
        # shifting everything right and causing spurious "bad x,y,z" errors.
        # pdbqt_normalize / canonicalize_receptor will clamp this automatically.
        if len(line) >= 66:
            bfac_raw = line[60:66].strip()
            try:
                if bfac_raw and float(bfac_raw) > _BFACTOR_WARN_LIMIT:
                    report.warnings.append(LintIssue(
                        i, "b_factor", bfac_raw, "BFAC_OVERFLOW",
                        f"B-factor {bfac_raw} > 99.99 fills the 6-char field with no leading "
                        "space; autogrid4 may misparse adjacent occupancy+bfac as one token. "
                        "Fix with: ultidock pdbqt canonicalize-receptor <file> -o <file>",
                    ))
            except ValueError:
                pass  # already caught by _check_numeric_field

        # ── AD atom type (must not be blank) ──────────────────────────────────
        atom_type = line[77:79].strip() if len(line) >= 79 else ""
        if not atom_type:
            report.errors.append(LintIssue(
                i, "atom_type", "",
                "MISSING", "AutoDock atom type (cols 78-79) is blank",
            ))

        # ── Line length (warning) ─────────────────────────────────────────────
        raw_len = len(line.rstrip("\r\n"))
        if raw_len < 79:
            report.warnings.append(LintIssue(
                i, "line", line,
                "SHORT_LINE", f"line length {raw_len} < 79 (minimum for atom type field)",
            ))

    return report


def _reformat_atom_line(line: str) -> str:
    """
    Rewrite all numeric columns in one ATOM/HETATM line using fixedfmt.
    Non-numeric fields (record, serial, names, chain, etc.) are passed through
    byte-for-byte.

    The line is padded to 80 characters so every slice is always valid.
    The atom type is always emitted as exactly 2 characters (right-justified,
    space-padded) so the output is always >= 79 chars and a subsequent
    pdbqt_check() call sees the field correctly.
    """
    line = line.ljust(80)

    x          = float(line[30:38])
    y          = float(line[38:46])
    z          = float(line[46:54])
    occupancy  = float(line[54:60])
    b_factor   = float(line[60:66])
    charge     = float(line[70:76])
    # Atom type: right-justify in 2 chars so "C" becomes "C " not just "C"
    ad_type    = line[77:79].strip()
    ad_type_2  = f"{ad_type:>2s}"   # e.g. "C " or "OA"

    return (
        line[0:30]               +  # record + serial + atom_name + alt_loc + res + chain + res_seq + icode
        fmt_coord(x)             +  # cols 31-38
        fmt_coord(y)             +  # cols 39-46
        fmt_coord(z)             +  # cols 47-54
        fmt_occupancy(occupancy) +  # cols 55-60
        fmt_bfactor(b_factor)    +  # cols 61-66
        line[66:70]              +  # 4-char gap
        fmt_charge(charge)       +  # cols 71-76
        " "                      +  # col 77
        ad_type_2                   # cols 78-79, always 2 chars
    )


def pdbqt_normalize(path: Path, out_path: Path) -> None:
    """
    Reformat ATOM/HETATM numeric columns in a ligand PDBQT using fixedfmt.

    v1 guarantees:
      - Serial numbers preserved
      - Atom order preserved (no reordering)
      - ROOT / BRANCH / ENDBRANCH / TORSDOF lines passed through verbatim
      - Line count preserved (splitlines -> join one-to-one)
      - LF-only line endings on output
      - ASCII encoding enforced

    Raises LintError if pdbqt_check() finds any errors.
    Raises FixedFmtError (propagated) if a reformatted value overflows.
    """
    # Lint first so we never write a partially-reformatted file. A half-fixed
    # PDBQT is worse than the original because it looks clean but still breaks.
    report = pdbqt_check(path)
    if not report.ok:
        msg = (
            f"{len(report.errors)} error(s) in {path.name}; "
            "run 'ultidock pdbqt check <file>' for details"
        )
        raise LintError(msg)

    out_lines: list[str] = []
    for line in path.read_bytes().decode("ascii").splitlines(keepends=False):
        rec = line[:6].strip()
        if rec in _ATOM_RECORDS:
            line = _reformat_atom_line(line)
        out_lines.append(line)

    # \n only — no \r\n — guarantees byte identity across OS
    content = ("\n".join(out_lines) + "\n").encode("ascii")
    out_path.write_bytes(content)


# ─────────────────────────────────────────────────────────────────────────────
# Part C — Receptor canonicalization
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ReceptorAtom:
    """Parsed representation of one ATOM or HETATM record in a receptor PDBQT."""
    serial:       int
    atom_name:    str      # 4-char field, exact original spacing preserved
    alt_loc:      str      # single char
    res_name:     str      # 3-char
    chain_id:     str      # 1-char
    res_seq:      int      # residue sequence number (integer)
    i_code:       str      # insertion code, 1-char (usually ' ')
    x: float
    y: float
    z: float
    occupancy:    float
    b_factor:     float
    charge:       float
    ad_atom_type: str      # 2-char AD atom type, right-justified

    # Original record keyword ("ATOM  " or "HETATM") — preserved verbatim
    record: str = "ATOM  "


def _parse_receptor_atoms(path: Path) -> list[ReceptorAtom]:
    """Parse all ATOM/HETATM lines from a receptor PDBQT into ReceptorAtom dataclasses."""
    atoms: list[ReceptorAtom] = []
    for line in path.read_bytes().decode("ascii").splitlines():
        rec = line[:6]
        if rec.strip() not in _ATOM_RECORDS:
            continue
        line = line.ljust(80)
        atoms.append(ReceptorAtom(
            record    = rec,
            serial    = int(line[6:11].strip() or "0"),
            atom_name = line[12:16],
            alt_loc   = line[16],
            res_name  = line[17:20].strip(),
            chain_id  = line[21],
            res_seq   = int(line[22:26].strip() or "0"),
            i_code    = line[26],
            x         = float(line[30:38]),
            y         = float(line[38:46]),
            z         = float(line[46:54]),
            occupancy = float(line[54:60]),
            b_factor  = float(line[60:66]),
            charge    = float(line[70:76]),
            ad_atom_type = line[77:79].strip(),
        ))
    return atoms


def _filter_altloc(atoms: list[ReceptorAtom]) -> list[ReceptorAtom]:
    """
    If any atom has a non-space AltLoc, keep only the first (alphabetically
    lowest) altloc letter; discard all others.

    Atoms with altloc=' ' (space) are always kept.
    Rationale: AutoDock-GPU does not model alternate conformations; both copies
    would confuse the grid generator and inflate atom counts.
    """
    present = sorted({a.alt_loc for a in atoms if a.alt_loc.strip()})
    if not present:
        return atoms
    keep_loc = present[0]   # e.g. 'A'
    result = [a for a in atoms if not a.alt_loc.strip() or a.alt_loc == keep_loc]
    # Normalise the kept altloc to space so output is clean
    for a in result:
        if a.alt_loc == keep_loc:
            a.alt_loc = " "
    return result


def _atom_sort_key(atom: ReceptorAtom) -> tuple:
    """
    Canonical sort key — pure chemistry, no dependence on original serials.

    Priority (ascending):
      1. chain_id       (str)  — alphabetical chain order
      2. res_seq        (int)  — residue sequence number
      3. res_name       (str)  — stable for HETATM ligands / waters at same res_seq
      4. atom_name      (str)  — N before CA before C before O (alphabetical)
      5. ad_atom_type   (str)  — final tiebreak
    """
    return (
        atom.chain_id,
        atom.res_seq,
        atom.res_name,
        atom.atom_name.strip(),
        atom.ad_atom_type,
    )


def _format_receptor_atom(atom: ReceptorAtom) -> str:
    """Emit one ATOM/HETATM line with all numeric fields through fixedfmt."""
    serial_str    = f"{atom.serial:5d}"
    atom_name_str = f"{atom.atom_name:<4s}"      # preserve original 4-char name field
    res_seq_str   = f"{atom.res_seq:4d}"

    return (
        f"{atom.record:<6s}"            # cols  1- 6
        f"{serial_str}"                 # cols  7-11
        f" "                            # col  12
        f"{atom_name_str}"              # cols 13-16
        f"{atom.alt_loc}"               # col  17
        f"{atom.res_name:<3s}"          # cols 18-20
        f" "                            # col  21
        f"{atom.chain_id}"              # col  22
        f"{res_seq_str}"                # cols 23-26
        f"{atom.i_code}"                # col  27
        f"   "                          # cols 28-30 (standard PDB gap)
        f"{fmt_coord(atom.x)}"          # cols 31-38
        f"{fmt_coord(atom.y)}"          # cols 39-46
        f"{fmt_coord(atom.z)}"          # cols 47-54
        f"{fmt_occupancy(atom.occupancy)}"  # cols 55-60
        f"{fmt_bfactor(atom.b_factor)}"     # cols 61-66
        f"    "                         # cols 67-70 (gap)
        f"{fmt_charge(atom.charge)}"    # cols 71-76
        f" "                            # col  77
        f"{atom.ad_atom_type:>2s}"      # cols 78-79
    ).rstrip()


def canonicalize_receptor(
    infile: Path,
    outfile: Path,
    *,
    timestamp: str | None = None,
) -> str:
    """
    Deterministically canonicalize a receptor PDBQT.

    Steps:
      1. Lint with pdbqt_check() — raises LintError if errors present.
      2. Parse all ATOM/HETATM lines into ReceptorAtom dataclasses.
      3. Filter altloc — keep only first letter.
      4. Sort by canonical chemistry key.
      5. Renumber serials 1, 2, 3 … contiguously.
      6. Emit with all numeric fields rewritten via fixedfmt.
      7. Write LF-only, ASCII-encoded output.
      8. Return SHA-256 hex digest of the output file.

    Parameters
    ----------
    infile    : path to the input receptor PDBQT
    outfile   : path to write the canonicalized output
    timestamp : if None, a UTC timestamp is embedded in the REMARK line;
                pass a fixed string (e.g. "FIXED") in tests to get a
                deterministic REMARK and thus a deterministic output hash.

    Returns
    -------
    SHA-256 hex digest (str) of the output file bytes.

    Raises
    ------
    LintError if pdbqt_check() reports any errors.
    """
    from molguard import __version__

    # ── 1. Lint ───────────────────────────────────────────────────────────────
    report = pdbqt_check(infile)
    if not report.ok:
        raise LintError(
            f"{len(report.errors)} error(s) in {infile.name}; "
            "run 'ultidock pdbqt check <file>' for details"
        )

    # ── 2. Parse ──────────────────────────────────────────────────────────────
    atoms = _parse_receptor_atoms(infile)

    # ── 3. AltLoc filter ─────────────────────────────────────────────────────
    n_before = len(atoms)
    atoms = _filter_altloc(atoms)
    n_discarded = n_before - len(atoms)
    if n_discarded:
        import warnings
        warnings.warn(
            f"{infile.name}: discarded {n_discarded} atom(s) with non-primary altloc",
            stacklevel=2,
        )

    # ── 4. Deterministic sort ─────────────────────────────────────────────────
    # Deterministic sort: pure chemistry key, no dependence on input serial order.
    # This is what makes the output byte-identical regardless of how the upstream
    # tool chose to number atoms.
    atoms.sort(key=_atom_sort_key)

    # ── 5. Renumber serials ───────────────────────────────────────────────────
    for new_serial, atom in enumerate(atoms, start=1):
        atom.serial = new_serial

    # ── 6 & 7. Emit ───────────────────────────────────────────────────────────
    ts = timestamp if timestamp is not None else _utc_now()
    header = f"REMARK Canonicalized by ultidock {__version__} on {ts}"
    lines = [header] + [_format_receptor_atom(a) for a in atoms]
    content = ("\n".join(lines) + "\n").encode("ascii")
    outfile.write_bytes(content)

    # ── 8. Return digest ──────────────────────────────────────────────────────
    return hashlib.sha256(content).hexdigest()


def _utc_now() -> str:
    """Return current UTC time as an ISO-8601 string (no microseconds)."""
    from datetime import datetime, timezone
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
