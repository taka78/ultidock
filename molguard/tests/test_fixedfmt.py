"""
tests/test_fixedfmt.py
~~~~~~~~~~~~~~~~~~~~~~
Unit tests for molguard.io.fixedfmt.

Coverage targets:
  - NaN / Inf rejection
  - Negative-zero normalisation
  - Overflow detection
  - No scientific notation in output
  - Exact width invariant
  - Decimal-point presence invariant
  - All named wrappers callable
"""

from __future__ import annotations

import math

import pytest

from molguard.io.fixedfmt import (
    FixedFmtError,
    fmt_bfactor,
    fmt_charge,
    fmt_coord,
    fmt_fixed,
    fmt_gpf_center,
    fmt_gpf_energy,
    fmt_gpf_spacing,
    fmt_occupancy,
)


# ── fmt_fixed core ────────────────────────────────────────────────────────────

class TestFmtFixed:

    def test_basic_positive(self):
        result = fmt_fixed(1.234, width=8, decimals=3)
        assert result == "   1.234"
        assert len(result) == 8

    def test_basic_negative(self):
        result = fmt_fixed(-1.234, width=8, decimals=3)
        assert result == "  -1.234"
        assert len(result) == 8

    def test_zero(self):
        result = fmt_fixed(0.0, width=8, decimals=3)
        assert result == "   0.000"

    def test_negative_zero_normalised(self):
        """IEEE-754 -0.0 must format identically to 0.0."""
        pos = fmt_fixed(0.0, width=8, decimals=3)
        neg = fmt_fixed(-0.0, width=8, decimals=3)
        assert pos == neg
        assert "-" not in neg

    def test_always_has_decimal(self):
        for v in [0.0, 1.0, -1.0, 999.9, 0.001]:
            assert "." in fmt_fixed(v, width=8, decimals=3)

    def test_never_exponent(self):
        """Values that a naive str(float) would render in scientific notation."""
        for v in [1e-10, 1e-5, -1e-7, 0.00001]:
            result = fmt_fixed(v, width=8, decimals=3)
            assert "e" not in result.lower(), f"exponent found in {result!r} for v={v}"

    def test_exact_width(self):
        for w in [6, 7, 8, 10]:
            result = fmt_fixed(1.5, width=w, decimals=2)
            assert len(result) == w

    def test_nan_raises(self):
        with pytest.raises(FixedFmtError, match="not finite"):
            fmt_fixed(float("nan"), width=8, decimals=3)

    def test_inf_raises(self):
        with pytest.raises(FixedFmtError, match="not finite"):
            fmt_fixed(float("inf"), width=8, decimals=3)

    def test_neg_inf_raises(self):
        with pytest.raises(FixedFmtError, match="not finite"):
            fmt_fixed(float("-inf"), width=8, decimals=3)

    def test_overflow_raises(self):
        """Value too large for the column must raise, not silently truncate."""
        # width=8, decimals=3 → max integer-part = 9999 (4 digits + sign)
        with pytest.raises(FixedFmtError, match="overflows"):
            fmt_fixed(99999.9, width=8, decimals=3)

    def test_boundary_value_ok(self):
        """Value just inside the overflow boundary should succeed."""
        result = fmt_fixed(999.999, width=8, decimals=3)
        assert "." in result
        assert len(result) == 8

    def test_small_positive_rounds_to_zero(self):
        """1e-10 rounds to 0.000 at 3-decimal precision — must not raise."""
        result = fmt_fixed(1e-10, width=8, decimals=3)
        assert result == "   0.000"


# ── Named wrappers ────────────────────────────────────────────────────────────

class TestNamedWrappers:

    def test_fmt_coord_width(self):
        assert len(fmt_coord(1.0)) == 8

    def test_fmt_charge_width(self):
        assert len(fmt_charge(0.123)) == 6

    def test_fmt_occupancy_width(self):
        assert len(fmt_occupancy(1.0)) == 6

    def test_fmt_bfactor_width(self):
        assert len(fmt_bfactor(0.0)) == 6

    def test_fmt_gpf_spacing_width(self):
        assert len(fmt_gpf_spacing(0.375)) == 7

    def test_fmt_gpf_center_width(self):
        assert len(fmt_gpf_center(10.5)) == 8

    def test_fmt_gpf_energy_width(self):
        assert len(fmt_gpf_energy(999.0)) == 8

    def test_fmt_gpf_spacing_example(self):
        """Standard AutoDock spacing 0.375 Å should format correctly."""
        assert fmt_gpf_spacing(0.375) == " 0.3750"

    def test_fmt_coord_example(self):
        assert fmt_coord(-1.234) == "  -1.234"

    def test_fmt_charge_example(self):
        assert fmt_charge(-0.412) == "-0.412"

    @pytest.mark.parametrize("wrapper,value", [
        (fmt_coord,       float("nan")),
        (fmt_charge,      float("inf")),
        (fmt_gpf_spacing, float("-inf")),
    ])
    def test_wrappers_reject_nonfinite(self, wrapper, value):
        with pytest.raises(FixedFmtError):
            wrapper(value)
