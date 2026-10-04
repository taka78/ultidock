"""Regression checks across the docking runner, pose handoff and MD entry point."""

import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from click.testing import CliRunner
import pytest

from cli import ultidock
from docking import run as pipeline
from docking.pose_provenance import record_docking_run, record_pose_provenance
from poses import select_poses
import workflow
from test_workflow import protocol_inputs


@pytest.fixture
def integrated(tmp_path, monkeypatch):
    protocol, directory, config = protocol_inputs(tmp_path)
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    monkeypatch.setattr(workflow, "check_dependencies", lambda *a, **kw: True)
    ligands = tmp_path / "ligands"
    ligands.mkdir()
    old_pose = next(directory.glob("*.pdbqt"))
    coordinates = "\n".join(line for line in old_pose.read_text().splitlines()
                            if line.startswith("ATOM")) + "\n"
    ligand = ligands / "ethanol.pdbqt"
    ligand.write_text(coordinates)
    # Older, better-scoring results must never win this invocation's ranking.
    old_pose.write_text(old_pose.read_text().replace("-5 0 0", "-90 0 0"))
    receptors = tmp_path / "staged_receptors"
    receptors.mkdir()
    receptor = receptors / "protein.pdbqt"
    receptor.write_text(Path(config["docking_receptor_pdbqt"]).read_text() + "REMARK staged\n")
    actual = SimpleNamespace(DOCKING_DIR=str(directory), LIGANDS_DIR=str(ligands),
                             MACRO_MOL_DIR=str(receptors), GPU_TYPE="CPU")
    events, jobs, commands = [], [], []

    def setup(*args):
        events.append("setup")
        return actual

    class Docking:
        def run(self):
            events.append("docking")
            fresh = directory / "protein__S2__ethanol_87654321.pdbqt"
            fresh.write_text("MODEL 3\nREMARK VINA RESULT: -6 0 0\n" + coordinates + "ENDMDL\n")
            record_pose_provenance(fresh, ligand, receptor, "protein", "S2", "vina")
            return [str(fresh)]

    def simulate(jobdir, through, reviewed, **kwargs):
        events.append("md")
        jobs.append((jobdir, through, reviewed))

    def dispatch(command, **kwargs):
        # Exercise the actual MD CLI parser, validation and snapshot creation.
        commands.append(command)
        assert Path(command[2]).name == "workflow.py"
        code = workflow.main(command[3:])
        if code and kwargs.get("check"):
            raise subprocess.CalledProcessError(code, command)
        return SimpleNamespace(returncode=code)

    monkeypatch.setattr(pipeline, "_ensure_config", setup)
    monkeypatch.setitem(sys.modules, "extract", SimpleNamespace(
        main=lambda **kw: events.append("extract")))
    monkeypatch.setitem(sys.modules, "dock_v02", SimpleNamespace(DockingProcessor=Docking))
    monkeypatch.setitem(sys.modules, "analyse_docking_results", SimpleNamespace(
        __file__="analysis.py", main=lambda: events.append("analysis")))
    monkeypatch.setattr(workflow, "run", simulate)
    monkeypatch.setattr(pipeline.subprocess, "run", dispatch)
    return SimpleNamespace(protocol=protocol, config=config, actual=actual, directory=directory,
                           root=root, events=events, jobs=jobs, commands=commands)


def test_pipeline_continues_with_only_current_poses_and_actual_staged_receptor(integrated):
    case = integrated
    pipeline.main(["--md-config", str(case.protocol)])
    assert case.events == ["setup", "extract", "docking", "analysis", "md"]
    jobdir, through, reviewed = case.jobs[0]
    assert through == "npt" and reviewed is False
    job = json.loads((jobdir / "job.json").read_text())
    pose = job["systems"][0]["pose"]
    assert (pose["score"], pose["model"], pose["binding_site"]) == (-6, 3, "2")
    assert (jobdir / job["config"]["docking_receptor_pdbqt"]).read_text().endswith("REMARK staged\n")
    handoff = json.loads(next(case.directory.glob("docking-run-*.md.json")).read_text())
    assert handoff["job_dir"] == str(jobdir) and handoff["status"] == "complete"
    assert (jobdir / "docking_manifest.json").exists()


def test_prepare_only_and_skipped_analysis_still_handoff(integrated):
    pipeline.main(["--md", str(integrated.protocol), "--md-through", "prepare", "--skip-analysis"])
    assert integrated.events == ["setup", "extract", "docking"]
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "prepared"
    assert (Path(status["job_dir"]) / "job.json").is_file()


def test_md_preflight_failure_preserves_docking_and_reports_skip(integrated, monkeypatch):
    monkeypatch.setattr(workflow, "check_dependencies", lambda *a, **kw: False)
    assert pipeline.main(["--md-config", str(integrated.protocol)]) == 0
    assert integrated.events == ["setup", "extract", "docking", "analysis"]
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "skipped" and "preflight" in status["reason"]


@pytest.mark.parametrize("fault", ["engine", "chemistry"])
def test_md_input_mismatch_does_not_discard_screening_results(integrated, fault):
    if fault == "engine":
        integrated.actual.GPU_TYPE = "CUDA"
    else:
        integrated.config["ligands"] = {}
        integrated.protocol.write_text(json.dumps(integrated.config))
    if fault == "engine":
        pipeline.main(["--md-config", str(integrated.protocol)])
    else:
        with pytest.raises(subprocess.CalledProcessError):
            pipeline.main(["--md-config", str(integrated.protocol)])
    assert integrated.events == ["setup", "extract", "docking", "analysis"]


def test_docking_failure_stops_analysis_and_md(integrated, monkeypatch):
    class BrokenDocking:
        def run(self):
            raise RuntimeError("Vina failed")
    monkeypatch.setitem(sys.modules, "dock_v02", SimpleNamespace(DockingProcessor=BrokenDocking))
    with pytest.raises(RuntimeError, match="Vina failed"):
        pipeline.main(["--md-config", str(integrated.protocol)])
    assert integrated.events == ["setup", "extract"]
    assert not list(integrated.directory.glob("docking-run-*.json"))


def test_docking_only_does_not_call_md(integrated):
    pipeline.main([])
    assert integrated.events == ["setup", "extract", "docking", "analysis"]
    assert not integrated.commands


def test_partial_screening_runs_md_for_successes_and_records_failures(integrated, monkeypatch):
    original = sys.modules["dock_v02"].DockingProcessor
    class PartialDocking(original):
        failures = [{"stage": "docking", "ligand_id": "failed_ligand",
                     "receptor_id": "protein", "binding_site": "1", "error": "Vina timed out"}]
    monkeypatch.setitem(sys.modules, "dock_v02", SimpleNamespace(DockingProcessor=PartialDocking))
    pipeline.main(["--md-config", str(integrated.protocol)])
    assert integrated.events[-1] == "md"
    manifest = next(path for path in integrated.directory.glob("docking-run-*.json")
                    if not path.name.endswith(".md.json"))
    status = json.loads(manifest.read_text())
    assert status["status"] == "partial" and status["failed_cases"] == 1
    assert "failed_ligand" in Path(status["failure_report"]).read_text()
    job = json.loads((integrated.jobs[0][0] / "job.json").read_text())
    assert [entry["pose"]["ligand"] for entry in job["systems"]] == ["ethanol"]


def test_all_failed_writes_report_and_never_uses_stale_results(integrated, monkeypatch):
    class FailedDocking:
        failures = [{"stage": "receptor-grids", "receptor_id": "protein", "error": "bad grids"}]
        def run(self):
            return []
    monkeypatch.setitem(sys.modules, "dock_v02", SimpleNamespace(DockingProcessor=FailedDocking))
    assert pipeline.main(["--md-config", str(integrated.protocol)]) == 1
    assert integrated.events == ["setup", "extract"]
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "skipped"
    manifest = next(path for path in integrated.directory.glob("docking-run-*.json")
                    if not path.name.endswith(".md.json"))
    report = json.loads(manifest.read_text())
    assert report["status"] == "failed" and report["outputs"] == []


def test_failed_md_receptor_does_not_borrow_other_receptor_results(integrated, monkeypatch):
    original = sys.modules["dock_v02"].DockingProcessor
    class OtherReceptor(original):
        failures = [{"stage": "receptor-grids", "receptor_id": "protein", "error": "bad grids"}]
        def run(self):
            outputs = super().run()
            sidecar = Path(outputs[0]).with_suffix(".md.json")
            record = json.loads(sidecar.read_text())
            record["receptor_id"] = "other"
            sidecar.write_text(json.dumps(record))
            return outputs
    monkeypatch.setitem(sys.modules, "dock_v02", SimpleNamespace(DockingProcessor=OtherReceptor))
    assert pipeline.main(["--md-config", str(integrated.protocol)]) == 0
    assert integrated.events[-1] == "analysis" and not integrated.jobs
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "skipped" and "requested MD receptor" in status["reason"]


def test_md_failure_keeps_job_path_and_returns_failure(integrated, monkeypatch):
    def failed_run(*args, **kwargs):
        raise RuntimeError("GROMACS failed")
    monkeypatch.setattr(workflow, "run", failed_run)
    with pytest.raises(subprocess.CalledProcessError):
        pipeline.main(["--md-config", str(integrated.protocol)])
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "failed" and "GROMACS failed" in status["error"]
    assert (Path(status["job_dir"]) / "job.json").is_file()


def test_final_md_preflight_failure_records_reason_and_preserves_docking(integrated):
    # The submitted receptor passes the early check, but the actual staged
    # receptor has a different coordinate frame when the final handoff checks it.
    receptor = Path(integrated.actual.MACRO_MOL_DIR) / "protein.pdbqt"
    receptor.write_text(receptor.read_text().replace("   1.000", "  21.000"))
    with pytest.raises(subprocess.CalledProcessError):
        pipeline.main(["--md-config", str(integrated.protocol)])
    assert integrated.events == ["setup", "extract", "docking", "analysis"]
    status = json.loads(next(integrated.directory.glob("docking-run-*.md.json")).read_text())
    assert status["status"] == "failed" and "coordinate frame" in status["error"]
    assert "job_dir" not in status and not integrated.jobs


def test_existing_results_run_prepares_and_executes(integrated):
    code = workflow.main(["run", "--config", str(integrated.protocol),
                          "--docking-dir", str(integrated.directory)])
    assert code == 0 and integrated.jobs[0][1] == "npt"


def test_new_jobs_cannot_claim_equilibration_review(integrated):
    code = workflow.main(["run", "--config", str(integrated.protocol),
                          "--docking-dir", str(integrated.directory), "--through", "production",
                          "--equilibration-reviewed"])
    assert code == 1 and not integrated.jobs


def test_manifest_detects_changed_outputs(tmp_path):
    _, directory, _ = protocol_inputs(tmp_path)
    pose = next(directory.glob("*.pdbqt"))
    manifest = record_docking_run(directory, [pose])
    pose.write_text(pose.read_text().replace("-5 0 0", "-8 0 0"))
    with pytest.raises(ValueError, match="changed or missing"):
        select_poses(directory, "protein", "vina", manifest=manifest)


@pytest.mark.parametrize("mode, options", [
    ("run", []), ("known-site", ["--center", "1,2,3"]),
    ("cavity", []), ("blind", []),
    ("p2rank", ["--receptor", "protein.pdbqt"]),
    ("fpocket", ["--receptor", "protein.pdbqt"]),
    ("run", ["p2rank", "--receptor", "protein.pdbqt"]),
])
def test_public_modes_keep_caller_relative_md_paths(tmp_path, monkeypatch, mode, options):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "protocol.json").write_text("{}")
    (tmp_path / "protein.pdbqt").write_text("fixture")
    seen = []
    monkeypatch.setattr(ultidock, "_docking_dir", lambda: tmp_path)
    monkeypatch.setattr(ultidock, "_run_python", lambda script, args, **kw: seen.append(args))
    monkeypatch.setattr(ultidock, "_run_pipeline_mode", lambda **kw: seen.append(kw["extra_args"]))
    monkeypatch.setattr(ultidock, "_create_pocket_sites", lambda *a: None)
    result = CliRunner().invoke(ultidock.cli, [mode, *options, "--md-config", "protocol.json",
        "--md-work-dir", "md-simulation/workspace/test", "--md-gmx", "./tools/gmx",
        "--ligands-dir=ligands"])
    assert result.exit_code == 0, result.output + str(result.exception)
    forwarded = seen[0]
    assert forwarded[forwarded.index("--md-config") + 1] == str(tmp_path / "protocol.json")
    assert forwarded[forwarded.index("--md-gmx") + 1] == str(tmp_path / "tools/gmx")
    assert f"--ligands-dir={tmp_path / 'ligands'}" in forwarded


def test_md_flags_are_visible_and_require_protocol():
    result = CliRunner().invoke(ultidock.cli, ["run", "--help"])
    assert result.exit_code == 0 and "--md-config" in result.output
    with pytest.raises(SystemExit) as exc:
        pipeline.main(["--md-through", "npt"])
    assert exc.value.code == 2
