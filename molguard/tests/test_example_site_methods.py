"""Example site selection must reach both pocket prediction and docking."""

from __future__ import annotations

import runpy
import subprocess
import sys
from unittest.mock import Mock

import pytest
from click.testing import CliRunner

from cli import ultidock
from examples import common


def test_example_stages_only_its_own_inputs(monkeypatch, tmp_path):
    monkeypatch.setattr(common, "DOCKING_DIR", tmp_path / "shared_tools")
    old_ligands = common.DOCKING_DIR / "LIGANDS_DIR"
    old_ligands.mkdir(parents=True)
    (old_ligands / "unrelated.pdbqt").write_text("old ligand", encoding="utf-8")
    example = tmp_path / "example"
    example.mkdir()
    receptor = example / "receptor.pdbqt"
    ligand = example / "ligand.pdbqt"
    receptor.write_text("receptor", encoding="utf-8")
    ligand.write_text("ligand", encoding="utf-8")

    first = common.stage_inputs(example, receptor, [ligand])
    second = common.stage_inputs(example, receptor, [ligand])

    assert first["workspace"] != second["workspace"]
    assert [path.name for path in first["ligands"].glob("*.pdbqt")] == ["ligand.pdbqt"]
    assert first["receptor"] == first["macro"] / "receptor.pdbqt"
    assert (old_ligands / "unrelated.pdbqt").exists()


@pytest.mark.parametrize("method", ["p2rank", "fpocket"])
def test_example_pocket_method_reaches_pipeline(monkeypatch, tmp_path, method):
    paths = {
        key: tmp_path / key
        for key in ("receptor", "ligands", "docking", "analysis", "vina", "autodock", "macro", "results")
    }
    pocket_boxes = Mock()
    pipeline = Mock()
    monkeypatch.setattr("docking.pocket_boxes.create_pocket_boxes", pocket_boxes)
    monkeypatch.setattr(common.subprocess, "run", pipeline)

    common.run_pipeline(paths, mode="gpu", site_method=method)

    assert pocket_boxes.call_args.kwargs["method"] == method
    assert pocket_boxes.call_args.kwargs["receptor_pdbqt"] == paths["receptor"]
    centers = paths["results"] / f"{method}-sites.tsv"
    assert pocket_boxes.call_args.kwargs["output_tsv"] == centers
    command = pipeline.call_args.args[0]
    assert command[command.index("--centers-tsv") + 1] == str(centers.resolve())
    assert command[command.index("--grid-mode") + 1] == "centers"
    assert "--skip-wget" in command


def test_default_example_sites_stay_in_workspace(monkeypatch, tmp_path):
    paths = {key: tmp_path / key for key in ("ligands", "docking", "analysis", "vina", "autodock", "macro", "results")}
    pipeline = Mock()
    monkeypatch.setattr(common.subprocess, "run", pipeline)

    common.run_pipeline(paths, mode="cpu")

    command = pipeline.call_args.args[0]
    assert command[command.index("--centers-tsv") + 1] == str((paths["results"] / "cav-emps-sites.tsv").resolve())
    assert command[command.index("--mode") + 1] == "cpu"


@pytest.mark.parametrize(
    ("extra_args", "expected_mode"),
    [(["p2rank"], "gpu"), (["p2rank", "--mode", "cpu"], "cpu")],
)
def test_sert_example_forwards_method(monkeypatch, tmp_path, extra_args, expected_mode):
    pipeline = Mock()
    monkeypatch.setitem(sys.modules, "common", common)
    monkeypatch.setattr(common, "stage_inputs", lambda *_: {"receptor": tmp_path / "receptor.pdbqt"})
    monkeypatch.setattr(common, "run_pipeline", pipeline)
    monkeypatch.setattr(sys, "argv", ["example-run.py", *extra_args])

    runpy.run_path(str(common.REPO_ROOT / "examples/sert-escitalopram/example-run.py"), run_name="__main__")

    assert pipeline.call_args.kwargs == {"mode": expected_mode, "site_method": "p2rank"}


@pytest.mark.parametrize(
    "extra_args",
    [["p2rank"], ["p2rank", "--mode", "cpu"]],
)
def test_example_cli_passes_method_to_runner(monkeypatch, extra_args):
    runner = Mock()
    monkeypatch.setattr(ultidock, "_run_python", runner)

    result = CliRunner().invoke(ultidock.cli, ["example", "run", "sert-escitalopram", *extra_args])

    assert result.exit_code == 0
    assert runner.call_args.args[1] == tuple(extra_args)


def test_pocket_module_imports_after_example_common():
    command = [
        sys.executable,
        "-c",
        "import sys; sys.path.insert(0, 'examples'); "
        "import common; import docking.pocket_boxes; "
        "assert common.__file__.endswith('examples/common.py')",
    ]
    subprocess.run(command, cwd=common.REPO_ROOT, check=True)


@pytest.mark.parametrize("method", ["cav-emps", "fpocket", "p2rank"])
def test_d2_stages_manifest_ligands_and_forwards_method(monkeypatch, tmp_path, method):
    paths = {"results": tmp_path / "results", "docking": tmp_path / "docking"}
    stage = Mock(return_value=paths)
    pipeline = Mock()
    monkeypatch.setitem(sys.modules, "common", common)
    monkeypatch.setattr(common, "stage_inputs", stage)
    monkeypatch.setattr(common, "run_pipeline", pipeline)
    monkeypatch.setattr(sys, "argv", ["example-run.py", method, "--mode", "cpu"])

    runpy.run_path(str(common.REPO_ROOT / "examples/d2-antipsychotics/example-run.py"),
                   run_name="__main__")

    root, receptor, ligands = stage.call_args.args
    assert root.name == "d2-antipsychotics"
    assert receptor.name == "6CM4-edited.pdbqt"
    assert [path.name for path in ligands] == [
        "haloperidol.pdbqt", "escitalopram-e.pdbqt", "morphine-e.pdbqt",
    ]
    pipeline.assert_called_once_with(paths, mode="cpu", site_method=method)


def test_d2_preview_does_not_stage_inputs(monkeypatch):
    stage = Mock()
    pipeline = Mock()
    monkeypatch.setitem(sys.modules, "common", common)
    monkeypatch.setattr(common, "stage_inputs", stage)
    monkeypatch.setattr(common, "run_pipeline", pipeline)
    monkeypatch.setattr(sys, "argv", ["example-run.py", "--dry-run"])

    runpy.run_path(str(common.REPO_ROOT / "examples/d2-antipsychotics/example-run.py"),
                   run_name="__main__")

    stage.assert_not_called()
    pipeline.assert_not_called()
