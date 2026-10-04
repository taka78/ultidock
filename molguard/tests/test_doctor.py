"""The main doctor must report docking and MD, even without an MD toolchain."""

from pathlib import Path
import shutil

from click.testing import CliRunner

from cli import ultidock


def test_doctor_reports_all_workflows_with_missing_md_tools(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[2]
    md = tmp_path / "md-simulation"
    md.mkdir()
    shutil.copy2(repo / "md-simulation/diagnostics.py", md / "diagnostics.py")
    vina = tmp_path / "docking/VINA_DIR/bin/vina"
    vina.parent.mkdir(parents=True)
    vina.write_text("#!/bin/sh\nexit 0\n")
    vina.chmod(0o755)
    monkeypatch.setattr(ultidock, "_repo_root", lambda: tmp_path)
    monkeypatch.setenv("PATH", "")
    result = CliRunner().invoke(ultidock.cli, ["doctor"])
    assert result.exit_code == 0, result.output
    for label in ("ultidock environment", "molguard", "python executable", "AutoDock Vina",
                  "receptor preparation", "Meeko receptor prep", "molecular dynamics", "GROMACS",
                  "ACPYPE", "Open Babel", "antechamber", "parmchk2", "tleap", "sqm",
                  "Python rdkit", "Python numpy", "Python scipy"):
        assert label in result.output
    assert str(vina.resolve()) in result.output
    assert "gmx not found; required for MD" in result.output
    assert f"MD workspace: {md / 'workspace'}" in result.output


def test_doctor_checks_custom_md_executable(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[2]
    md = tmp_path / "md-simulation"
    md.mkdir()
    shutil.copy2(repo / "md-simulation/diagnostics.py", md / "diagnostics.py")
    gmx = tmp_path / "custom-gromacs"
    gmx.write_text("#!/bin/sh\nexit 99\n")
    gmx.chmod(0o755)
    monkeypatch.setattr(ultidock, "_repo_root", lambda: tmp_path)
    monkeypatch.setenv("PATH", "")
    result = CliRunner().invoke(ultidock.cli, ["doctor", "--gmx", str(gmx)])
    assert result.exit_code == 0, result.output
    assert any(line.split() == ["[OK]", "GROMACS", str(gmx)] for line in result.output.splitlines())
    assert "gmx not found" not in result.output
