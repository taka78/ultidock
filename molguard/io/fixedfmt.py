"""
ultidock.io.fixedfmt
~~~~~~~~~~~~~~~~~~~~
Fixed-width, locale-safe float formatting for AutoDock-era column parsers.

All AutoGrid / AutoDock4 file formats use Fortran-style column ranges.
A single misplaced character (exponent, missing decimal, extra space) silently
corrupts the parse.  Every float that touches those files MUST go through one
of the functions in this module.

Design invariants (enforced on every call):
  - Never emits scientific notation (no 'e' / 'E')
  - Always contains exactly one decimal point
  - Output width is always exactly `width` characters
  - Decimal separator is always '.' regardless of OS locale
  - NaN and Inf are rejected before any formatting attempt
  - Negative zero (-0.0) is normalised to positive zero (0.0)
"""

from __future__ import annotations

import math

__all__ = [
    "FixedFmtError",
    "fmt_fixed",
    "fmt_coord",
    "fmt_charge",
    "fmt_occupancy",
    "fmt_bfactor",
    "fmt_gpf_spacing",
    "fmt_gpf_center",
    "fmt_gpf_energy",
]


class FixedFmtError(ValueError):
    """Raised when a value cannot be represented in the target fixed-width format."""


def fmt_fixed(
    value: float,
    *,
    width: int,
    decimals: int,
    field_name: str = "value",
) -> str:
    """
    Core formatter.  Returns a right-justified, fixed-width decimal string.

    Parameters
    ----------
    value      : float to format
    width      : total character width of the returned string (including sign and '.')
    decimals   : number of digits after the decimal point
    field_name : label used in error messages to identify which field failed

    Returns
    -------
    str of exactly `width` characters, right-justified, no exponent notation.

    Raises
    ------
    FixedFmtError
        If value is NaN, Inf, or its absolute value overflows the integer-part
        of the field (i.e., it does not fit in width chars with `decimals` decimals).
    """
    # ── 1. Reject non-finite ──────────────────────────────────────────────────
    # Python float() can silently produce nan/inf from certain inputs; we catch
    # them here before they reach the file and confuse Fortran parsers downstream.
    if not math.isfinite(value):
        raise FixedFmtError(
            f"{field_name}={value!r} is not finite — refusing to write to file"
        )

    # ── 2. Normalise negative zero (IEEE 754: -0.0 == 0.0 but formats differ) ─
    # f"{-0.0:8.3f}" produces "  -0.000" which is technically wrong for a charge
    # of zero and triggers spurious diff noise between runs.
    if value == 0.0:
        value = 0.0  # re-assign to strip the sign bit cleanly

    # ── 3. Overflow guard ─────────────────────────────────────────────────────
    # AutoGrid uses a fixed field width; if the integer part overflows it, the
    # column boundary shifts and everything to the right is misread silently.
    # We raise here so the caller knows exactly which field and value caused it.
    max_abs = 10 ** (width - decimals - 1)
    if abs(value) >= max_abs:
        raise FixedFmtError(
            f"{field_name}={value} overflows width={width} with decimals={decimals} "
            f"(max absolute value is {max_abs - 10**-decimals:.{decimals}f})"
        )

    # ── 4. Format — Python's f-string is always locale-independent ────────────
    result = f"{value:{width}.{decimals}f}"

    # ── 5. Self-check (cheap; catches any future Python regressions) ──────────
    assert len(result) == width, f"width mismatch: got {len(result)}, want {width}: {result!r}"
    assert "." in result, f"no decimal point in result: {result!r}"
    assert "e" not in result and "E" not in result, f"exponent leaked into result: {result!r}"

    return result


# ── Convenience wrappers ──────────────────────────────────────────────────────
# Width/precision pinned to the AutoDock PDBQT / GPF column specification.
# Never call fmt_fixed() ad-hoc in production code; use these named wrappers.

def fmt_coord(value: float) -> str:
    """ATOM/HETATM X, Y, Z coordinate column: 8 chars, 3 decimals.
    Example output: '  -1.234'  or '  99.999'"""
    return fmt_fixed(value, width=8, decimals=3, field_name="coord")


def fmt_charge(value: float) -> str:
    """PDBQT partial-charge column (cols 70-75): 6 chars, 3 decimals.
    Example output: '-0.412'  or ' 0.523'"""
    return fmt_fixed(value, width=6, decimals=3, field_name="charge")


def fmt_occupancy(value: float) -> str:
    """ATOM occupancy column (cols 55-60): 6 chars, 2 decimals.
    Example output: '  1.00'"""
    return fmt_fixed(value, width=6, decimals=2, field_name="occupancy")


def fmt_bfactor(value: float) -> str:
    """ATOM B-factor column (cols 61-66): 6 chars, 2 decimals.
    Example output: '  0.00'"""
    return fmt_fixed(value, width=6, decimals=2, field_name="bfactor")


def fmt_gpf_spacing(value: float) -> str:
    """GPF 'spacing' directive: 7 chars, 4 decimals.
    Example output: ' 0.3750'"""
    return fmt_fixed(value, width=7, decimals=4, field_name="spacing")


def fmt_gpf_center(value: float) -> str:
    """GPF 'gridcenter' coordinate: 8 chars, 3 decimals (same as PDBQT coord)."""
    return fmt_fixed(value, width=8, decimals=3, field_name="center")


def fmt_gpf_energy(value: float) -> str:
    """GPF energy-like fields (smooth, dielectric cutoff): 8 chars, 3 decimals."""
    return fmt_fixed(value, width=8, decimals=3, field_name="energy")
