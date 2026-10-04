"""Example site selection must reach both pocket prediction and docking."""

from __future__ import annotations

import runpy
import subprocess
import sys
from pathlib import Path
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
    [(["p2rank"], "auto"), (["p2rank", "--mode", "cpu"], "cpu")],
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


@pytest.mark.parametrize("name", ["sert-escitalopram", "gabaa-benzos"])
def test_docking_examples_forward_md_options(monkeypatch, tmp_path, name):
    pipeline = Mock()
    stage = Mock(return_value={"receptor": tmp_path / "protein.pdbqt"})
    monkeypatch.chdir(tmp_path)
    monkeypatch.setitem(sys.modules, "common", common)
    monkeypatch.setattr(common, "stage_inputs", stage)
    monkeypatch.setattr(common, "run_pipeline", pipeline)
    monkeypatch.setattr(sys, "argv", ["example-run.py", "--md-config", "protocol.json",
        "--md-through", "prepare", "--md-work-dir", "md-simulation/workspace/job",
        "--output-dir", "screen", "--dry-run"])
    runpy.run_path(str(common.REPO_ROOT / f"examples/{name}/example-run.py"), run_name="__main__")
    args = pipeline.call_args.kwargs
    assert args["mode"] == "auto" and args["dry_run"]
    assert args["extra_args"][:2] == ["--md-config", str(tmp_path / "protocol.json")]
    assert "--md-through" in args["extra_args"]
    assert stage.call_args.kwargs["output_dir"] == Path("screen")


def test_example_cli_resolves_paths_before_changing_directory(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    runner = Mock()
    monkeypatch.setattr(ultidock, "_run_python", runner)
    result = CliRunner().invoke(ultidock.cli, ["example", "run", "gabaa-8dd2-cav-emps",
        "--pdb", "input/8DD2.pdb", "--output-dir=study"])
    assert result.exit_code == 0
    assert runner.call_args.args[1] == ("--pdb", str(tmp_path / "input/8DD2.pdb"),
                                        f"--output-dir={tmp_path / 'study'}")


def test_example_preview_does_not_run_or_install_pocket_tools(monkeypatch, tmp_path):
    paths = {key: tmp_path / key for key in
             ("receptor", "ligands", "docking", "analysis", "vina", "autodock", "macro", "results")}
    pocket = Mock(side_effect=AssertionError("preview launched a pocket tool"))
    command = Mock(side_effect=AssertionError("preview launched the pipeline"))
    monkeypatch.setattr("docking.pocket_boxes.create_pocket_boxes", pocket)
    monkeypatch.setattr(common.subprocess, "run", command)
    common.run_pipeline(paths, site_method="p2rank", dry_run=True)


def test_example_refuses_to_overwrite_existing_workspace(tmp_path):
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "results.csv"
    marker.write_text("previous results")
    with pytest.raises(FileExistsError):
        common._ensure_workspace(tmp_path, output)
    assert marker.read_text() == "previous results"
