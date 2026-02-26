"""
tests/test_pdbqt.py
~~~~~~~~~~~~~~~~~~~
Tests for ultidock.io.pdbqt — lint, normalize, and receptor canonicalization.

Test classes:
  TestPdbqtCheck               — pdbqt_check() regression tests
  TestPdbqtNormalize           — pdbqt_normalize() determinism & invariants
  TestReceptorCanonicalization — canonicalize_receptor() idempotency & invariants
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from molguard.io.pdbqt import (
    LintError,
    canonicalize_receptor,
    pdbqt_check,
    pdbqt_normalize,
)

FIXTURES = Path(__file__).parent / "fixtures"


# ── pdbqt_check ───────────────────────────────────────────────────────────────

class TestPdbqtCheck:

    def test_good_ligand_passes(self, good_ligand):
        report = pdbqt_check(good_ligand)
        assert report.ok, f"Unexpected errors: {report.errors}"

    def test_good_receptor_passes(self, good_receptor):
        report = pdbqt_check(good_receptor)
        assert report.ok, f"Unexpected errors: {report.errors}"

    # ── Regression: decimal point vanished ───────────────────────────────────

    def test_decimal_vanished_detected(self, bad_decimal):
        """
        Regression test for the 'decimal point vanished' bug.
        A coordinate like '  -1234' (no '.') must trigger a NO_DECIMAL error.
        """
        report = pdbqt_check(bad_decimal)
        assert not report.ok
        codes = {e.code for e in report.errors}
        assert "NO_DECIMAL" in codes, (
            f"Expected NO_DECIMAL error; got codes: {codes}\n"
            f"Errors: {report.errors}"
        )

    def test_decimal_vanished_on_correct_column(self, bad_decimal):
        """The NO_DECIMAL error must reference the X column (col 30-38)."""
        report = pdbqt_check(bad_decimal)
        x_errors = [e for e in report.errors if e.code == "NO_DECIMAL" and e.column == "X"]
        assert x_errors, "NO_DECIMAL not reported on the X column"

    # ── Regression: exponent notation ────────────────────────────────────────

    def test_exponent_in_charge_detected(self, exponent_field):
        """A charge field containing '1.23e-04' must trigger an EXPONENT error."""
        report = pdbqt_check(exponent_field)
        assert not report.ok
        codes = {e.code for e in report.errors}
        assert "EXPONENT" in codes, (
            f"Expected EXPONENT error; got codes: {codes}\n"
            f"Errors: {report.errors}"
        )

    def test_exponent_on_correct_column(self, exponent_field):
        report = pdbqt_check(exponent_field)
        exp_errors = [e for e in report.errors if e.code == "EXPONENT" and e.column == "charge"]
        assert exp_errors, "EXPONENT not reported on the charge column"

    # ── Non-ASCII ─────────────────────────────────────────────────────────────

    def test_non_ascii_detected(self, tmp_path):
        bad = tmp_path / "nonascii.pdbqt"
        bad.write_bytes(b"ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00   \xff0.50 C\n")
        report = pdbqt_check(bad)
        assert not report.ok
        assert any(e.code == "NON_ASCII" for e in report.errors)

    # ── Missing atom type ─────────────────────────────────────────────────────

    def test_missing_atom_type_detected(self, tmp_path):
        atom_line = "ATOM      1  C1  LIG A   1      -1.234   2.567  -0.001  1.00  0.00    -0.101\n"
        pdbqt = tmp_path / "no_type.pdbqt"
        pdbqt.write_text(atom_line, encoding="ascii")
        report = pdbqt_check(pdbqt)
        assert not report.ok
        assert any(e.code == "MISSING" and e.column == "atom_type" for e in report.errors)

    # ── Fail-slow: collect all errors ─────────────────────────────────────────

    def test_multiple_errors_collected(self, tmp_path):
        """All lines must be checked before returning; no early exit on first error."""
        lines = [
            "ATOM      1  C1  LIG A   1      -1234   2.567  -0.001  1.00  0.00    -0.101 C\n",  # NO_DECIMAL
            "ATOM      2  N1  LIG A   1       0.000   0.000   0.000  1.00  0.00  1.0e-3 N\n",   # EXPONENT
        ]
        pdbqt = tmp_path / "multi.pdbqt"
        pdbqt.write_text("".join(lines), encoding="ascii")
        report = pdbqt_check(pdbqt)
        codes = {e.code for e in report.errors}
        assert "NO_DECIMAL" in codes
        assert "EXPONENT" in codes


# ── pdbqt_normalize ───────────────────────────────────────────────────────────

class TestPdbqtNormalize:

    def test_normalize_succeeds_on_good_ligand(self, good_ligand, tmp_path):
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_normalize_raises_on_bad_input(self, bad_decimal, tmp_path):
        out = tmp_path / "out.pdbqt"
        with pytest.raises(LintError):
            pdbqt_normalize(bad_decimal, out)

    # ── Determinism ───────────────────────────────────────────────────────────

    def test_normalize_deterministic(self, good_ligand, tmp_path):
        """Same input → byte-identical output on two successive runs."""
        out_a = tmp_path / "a.pdbqt"
        out_b = tmp_path / "b.pdbqt"
        pdbqt_normalize(good_ligand, out_a)
        pdbqt_normalize(good_ligand, out_b)
        assert out_a.read_bytes() == out_b.read_bytes(), (
            "pdbqt_normalize produced different output on two runs from the same input"
        )

    def test_normalize_hash_stable(self, good_ligand, tmp_path):
        out_a = tmp_path / "a.pdbqt"
        out_b = tmp_path / "b.pdbqt"
        pdbqt_normalize(good_ligand, out_a)
        pdbqt_normalize(good_ligand, out_b)
        h_a = hashlib.sha256(out_a.read_bytes()).hexdigest()
        h_b = hashlib.sha256(out_b.read_bytes()).hexdigest()
        assert h_a == h_b

    # ── Structural invariants ─────────────────────────────────────────────────

    def test_torsion_tree_preserved(self, good_ligand, tmp_path):
        """ROOT / BRANCH / ENDBRANCH / TORSDOF lines must be unchanged."""
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)

        src_tree = [
            l for l in good_ligand.read_text("ascii").splitlines()
            if l[:6].strip() in {"ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"}
        ]
        dst_tree = [
            l for l in out.read_text("ascii").splitlines()
            if l[:6].strip() in {"ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"}
        ]
        assert src_tree == dst_tree, "Torsion tree lines were modified"

    def test_line_count_preserved(self, good_ligand, tmp_path):
        """Output must have the same number of lines as input."""
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)
        src_n = len(good_ligand.read_text("ascii").splitlines())
        dst_n = len(out.read_text("ascii").splitlines())
        assert src_n == dst_n

    def test_serial_numbers_preserved(self, good_ligand, tmp_path):
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)

        def serials(text: str) -> list[str]:
            return [
                line[6:11]
                for line in text.splitlines()
                if line[:6].strip() in {"ATOM", "HETATM"}
            ]

        src_serials = serials(good_ligand.read_text("ascii"))
        dst_serials = serials(out.read_text("ascii"))
        assert src_serials == dst_serials

    def test_lf_only_line_endings(self, good_ligand, tmp_path):
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)
        raw = out.read_bytes()
        assert b"\r\n" not in raw, "Output contains CRLF line endings"

    def test_output_is_ascii(self, good_ligand, tmp_path):
        out = tmp_path / "out.pdbqt"
        pdbqt_normalize(good_ligand, out)
        out.read_bytes().decode("ascii")   # raises if non-ASCII


# ── canonicalize_receptor ─────────────────────────────────────────────────────

class TestReceptorCanonicalization:

    def test_canonicalize_succeeds(self, good_receptor, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        digest = canonicalize_receptor(good_receptor, out, timestamp="FIXED")
        assert out.exists()
        assert isinstance(digest, str) and len(digest) == 64  # SHA-256 hex

    def test_raises_on_bad_receptor(self, bad_decimal, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        with pytest.raises(LintError):
            canonicalize_receptor(bad_decimal, out)

    # ── Idempotency ───────────────────────────────────────────────────────────

    def test_canonicalize_idempotent(self, good_receptor, tmp_path):
        """Running canonicalize_receptor twice must produce byte-identical output."""
        out1 = tmp_path / "r1.pdbqt"
        out2 = tmp_path / "r2.pdbqt"
        canonicalize_receptor(good_receptor, out1, timestamp="FIXED")
        canonicalize_receptor(out1, out2, timestamp="FIXED")
        assert out1.read_bytes() == out2.read_bytes(), (
            "canonicalize_receptor is not idempotent"
        )

    # ── Determinism ───────────────────────────────────────────────────────────

    def test_canonicalize_hash_stable(self, good_receptor, tmp_path):
        out_a = tmp_path / "a.pdbqt"
        out_b = tmp_path / "b.pdbqt"
        h_a = canonicalize_receptor(good_receptor, out_a, timestamp="FIXED")
        h_b = canonicalize_receptor(good_receptor, out_b, timestamp="FIXED")
        assert h_a == h_b

    # ── Structural invariants ─────────────────────────────────────────────────

    def test_output_is_ascii(self, good_receptor, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        canonicalize_receptor(good_receptor, out, timestamp="FIXED")
        out.read_bytes().decode("ascii")

    def test_lf_only_endings(self, good_receptor, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        canonicalize_receptor(good_receptor, out, timestamp="FIXED")
        assert b"\r\n" not in out.read_bytes()

    def test_serials_contiguous(self, good_receptor, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        canonicalize_receptor(good_receptor, out, timestamp="FIXED")
        serials = []
        for line in out.read_text("ascii").splitlines():
            rec = line[:6].strip()
            if rec in {"ATOM", "HETATM"}:
                serials.append(int(line[6:11]))
        assert serials == list(range(1, len(serials) + 1)), (
            f"Serials not contiguous: {serials}"
        )

    def test_remark_line_present(self, good_receptor, tmp_path):
        out = tmp_path / "canonical.pdbqt"
        canonicalize_receptor(good_receptor, out, timestamp="FIXED")
        first_line = out.read_text("ascii").splitlines()[0]
        assert first_line.startswith("REMARK Canonicalized by ultidock")

    def test_atom_count_preserved(self, good_receptor, tmp_path):
        """All atoms must be present in output (no drops, assuming no altloc)."""
        out = tmp_path / "canonical.pdbqt"
        canonicalize_receptor(good_receptor, out, timestamp="FIXED")

        def count_atoms(p: Path) -> int:
            return sum(
                1 for l in p.read_text("ascii").splitlines()
                if l[:6].strip() in {"ATOM", "HETATM"}
            )

        assert count_atoms(good_receptor) == count_atoms(out)
