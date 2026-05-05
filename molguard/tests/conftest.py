"""
tests/conftest.py — shared fixtures for the MolGuard test suite.
"""
from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def good_ligand() -> Path:
    return FIXTURES_DIR / "good_ligand.pdbqt"


@pytest.fixture
def good_receptor() -> Path:
    return FIXTURES_DIR / "good_receptor.pdbqt"


@pytest.fixture
def bad_decimal() -> Path:
    return FIXTURES_DIR / "bad_decimal.pdbqt"


@pytest.fixture
def exponent_field() -> Path:
    return FIXTURES_DIR / "exponent_field.pdbqt"
