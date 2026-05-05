"""
tests/test_grids_check.py
~~~~~~~~~~~~~~~~~~~~~~~~~
Tests for molguard.grids.check — .fld / .map sanity checker.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from molguard.grids.check import GridSanityError, check_fld


# ── Helpers to build fake .fld + .map files ───────────────────────────────────

def _write_fake_map(directory: Path, filename: str, values: list[float]) -> Path:
    """Write a minimal AutoGrid-style .map file."""
    path = directory / filename
    lines = [
        "GRID_PARAMETER_FILE fake.gpf",
        "GRID_DATA_FILE fake.fld",
        "MACROMOLECULE fake.pdbqt",
    ]
    lines += [f"{v:.6f}" for v in values]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def _write_fake_fld(
    directory: Path,
    atom_types: list[str],
    stem: str = "receptor",
) -> Path:
    """Write a minimal .fld file referencing <stem>.<type>.map for each type."""
    fld_path = directory / f"{stem}.maps.fld"
    lines = [
        "# AVS field file",
        f"dim1=10",
        f"dim2=10",
        f"dim3=10",
        f"#SPACING 0.375",
        f"#CENTER 0.000 0.000 0.000",
        f"#NELEMENTS 10 10 10",
    ]
    for idx, at in enumerate(atom_types, start=1):
        map_file = f"{stem}.{at}.map"
        lines.append(f"variable {idx} file={map_file} filetype=ascii skip=6")
        lines.append(f"label={at}-affinity")
    fld_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return fld_path


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestCheckFld:

    def test_all_ok(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C", "OA"])
        _write_fake_map(tmp_path, "receptor.C.map",  [1.5, -0.3, 2.1, 0.0, -1.0])
        _write_fake_map(tmp_path, "receptor.OA.map", [0.5,  2.3, -1.1])
        report = check_fld(fld)
        assert report.ok

    def test_all_zero_map_detected(self, tmp_path):
        """An all-zero .map should trigger a clear error."""
        fld = _write_fake_fld(tmp_path, ["C", "OA"])
        _write_fake_map(tmp_path, "receptor.C.map",  [1.5, -0.3, 2.1])
        _write_fake_map(tmp_path, "receptor.OA.map", [0.0] * 50)
        report = check_fld(fld)
        assert not report.ok
        oa = next(r for r in report.map_results if r.atom_type == "OA")
        assert not oa.ok
        assert any("0.0" in e or "zero" in e.lower() for e in oa.errors), \
            f"Expected all-zero error; got: {oa.errors}"

    def test_missing_map_detected(self, tmp_path):
        """A referenced .map file that does not exist must be flagged."""
        fld = _write_fake_fld(tmp_path, ["C"])
        # Deliberately do NOT write receptor.C.map
        report = check_fld(fld)
        assert not report.ok
        c = next(r for r in report.map_results if r.atom_type == "C")
        assert not c.ok
        assert any("not found" in e.lower() or "exist" in e.lower() for e in c.errors), \
            f"Expected 'not found' error; got: {c.errors}"

    def test_type_mismatch_detected(self, tmp_path):
        """If expected_types doesn't match actual types, report an error."""
        fld = _write_fake_fld(tmp_path, ["C", "OA"])
        _write_fake_map(tmp_path, "receptor.C.map",  [1.0])
        _write_fake_map(tmp_path, "receptor.OA.map", [1.0])
        # Pretend we expected C, HD, OA — HD is missing
        report = check_fld(fld, expected_types=["C", "HD", "OA"])
        assert not report.ok
        assert any("mismatch" in e.lower() or "missing" in e.lower()
                   for e in report.errors)

    def test_expected_types_match_passes(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C", "OA"])
        _write_fake_map(tmp_path, "receptor.C.map",  [1.5])
        _write_fake_map(tmp_path, "receptor.OA.map", [0.5])
        report = check_fld(fld, expected_types=["C", "OA"])
        assert report.ok

    def test_fail_fast_raises(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C"])
        # Missing map → fail
        with pytest.raises(GridSanityError):
            check_fld(fld).fail_fast()

    def test_ok_report_does_not_raise(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C"])
        _write_fake_map(tmp_path, "receptor.C.map", [1.0, -0.5, 2.3])
        check_fld(fld).fail_fast()   # must not raise

    def test_sampled_range_recorded(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C"])
        _write_fake_map(tmp_path, "receptor.C.map", [-5.0, 0.0, 10.5])
        report = check_fld(fld)
        c = report.map_results[0]
        assert c.sampled_min == pytest.approx(-5.0)
        assert c.sampled_max == pytest.approx(10.5)

    def test_absurd_range_warns(self, tmp_path):
        fld = _write_fake_fld(tmp_path, ["C"])
        _write_fake_map(tmp_path, "receptor.C.map", [200000.0, -200000.0])
        report = check_fld(fld)
        c = report.map_results[0]
        # Should be a warning, not an error
        assert any("range" in w.lower() or "extreme" in w.lower() or "absurd" in w.lower()
                   for w in c.warnings), f"Expected absurd-range warning; got: {c.warnings}"

    def test_corrupt_fld_reports_error(self, tmp_path):
        """A .fld with no variable entries should produce a file-level error."""
        fld = tmp_path / "bad.maps.fld"
        fld.write_text("# Nothing useful here\ndim1=10\ndim2=10\ndim3=10\n", encoding="ascii")
        report = check_fld(fld)
        assert not report.ok
        assert report.errors  # fld-level error expected
