import builtins
from pathlib import Path
import runpy


def test_diagnostics_can_report_missing_scientific_packages(monkeypatch, capsys):
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name.split(".")[0] in {"numpy", "scipy", "rdkit"}:
            raise AssertionError("Diagnostics must not import the scientific stack")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    diagnostics = runpy.run_path(str(Path(__file__).resolve().parents[1] / "diagnostics.py"))
    monkeypatch.setattr(diagnostics["shutil"], "which", lambda _: None)
    monkeypatch.setattr(diagnostics["importlib"].util, "find_spec", lambda _: None)
    assert not diagnostics["check_dependencies"]()
    output = capsys.readouterr().out
    assert "Python numpy" in output and "Python scipy" in output and "Python rdkit" in output
    assert "MD dependencies: incomplete" in output
