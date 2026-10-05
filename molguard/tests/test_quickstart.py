"""The teacher must defer real docking until the final Next or explicit --yes."""

import importlib.util
from pathlib import Path
import sys
from unittest.mock import Mock

import pytest


@pytest.fixture
def teacher(monkeypatch):
    path = Path(__file__).resolve().parents[2] / "examples/quickstart/example-run.py"
    spec = importlib.util.spec_from_file_location("quickstart_teacher", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "environment", lambda *_: (["Tools available"], []))
    monkeypatch.setattr(sys, "stdin", Mock(isatty=lambda: True))
    monkeypatch.setattr(sys, "argv", [str(path)])
    monkeypatch.setattr(module.subprocess, "run", Mock(return_value=Mock(returncode=0)))
    return module


def test_preview_never_runs_pipeline(teacher, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["teacher", "--dry-run"])
    monkeypatch.setattr("builtins.input", Mock(side_effect=AssertionError("must not prompt")))
    teacher.main()
    teacher.subprocess.run.assert_not_called()
    output = capsys.readouterr().out
    assert "Step 5 of 5" in output
    assert "d2-antipsychotics" in output
    assert "No files were staged" in output


def test_next_back_and_exit_before_launch(teacher, monkeypatch, capsys):
    answers = iter(["", "b", "", "", "", "", "q"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    teacher.main()
    teacher.subprocess.run.assert_not_called()
    output = capsys.readouterr().out
    assert output.count("Step 1 of 5") == 2
    assert "Step 5 of 5" in output


def test_final_next_launches_exact_example(teacher, monkeypatch):
    prompts = []

    def advance(prompt):
        teacher.subprocess.run.assert_not_called()
        prompts.append(prompt)
        return ""

    monkeypatch.setattr("builtins.input", advance)
    teacher.main()
    assert len(prompts) == 5
    teacher.subprocess.run.assert_called_once_with(
        [sys.executable, "-u", str(teacher.EXAMPLES / "d2-antipsychotics/example-run.py"),
         "cav-emps", "--mode", "cpu"],
        cwd=teacher.EXAMPLES / "d2-antipsychotics",
    )


def test_noninteractive_requires_explicit_run_flag(teacher, monkeypatch):
    monkeypatch.setattr(sys, "stdin", Mock(isatty=lambda: False))
    with pytest.raises(SystemExit) as exc:
        teacher.main()
    assert exc.value.code == 2
    teacher.subprocess.run.assert_not_called()


def test_yes_forwards_choices_and_failure(teacher, monkeypatch):
    monkeypatch.setattr(sys, "stdin", Mock(isatty=lambda: False))
    monkeypatch.setattr(sys, "argv", ["teacher", "--yes", "--example", "sert-escitalopram",
                                     "--mode", "gpu", "--site-method", "fpocket"])
    teacher.subprocess.run.return_value.returncode = 7
    with pytest.raises(SystemExit) as exc:
        teacher.main()
    assert exc.value.code == 7
    command = teacher.subprocess.run.call_args.args[0]
    assert command[-3:] == ["fpocket", "--mode", "gpu"]
    assert "sert-escitalopram" in command[2]


def test_missing_tools_block_unattended_launch(teacher, monkeypatch):
    monkeypatch.setattr(teacher, "environment", lambda *_: ([], ["Vina missing"]))
    monkeypatch.setattr(sys, "argv", ["teacher", "--yes"])
    with pytest.raises(SystemExit) as exc:
        teacher.main()
    assert exc.value.code == 1
    teacher.subprocess.run.assert_not_called()


@pytest.mark.parametrize("error", [EOFError, KeyboardInterrupt])
def test_interrupted_prompt_exits_without_launch(teacher, monkeypatch, error):
    monkeypatch.setattr("builtins.input", Mock(side_effect=error))
    teacher.main()
    teacher.subprocess.run.assert_not_called()
