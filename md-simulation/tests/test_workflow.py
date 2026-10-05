import json
import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from poses import digest
from protocol import load_protocol
import workflow
from test_handoff import atom, ethanol_source


def protocol_inputs(tmp_path):
    receptor = tmp_path / "protein.pdb"
    receptor.write_text(atom(1, name="CA") + "\n")
    docked = tmp_path / "protein.pdbqt"
    docked.write_text(receptor.read_text())
    sdf = tmp_path / "ligand.sdf"
    ethanol = ethanol_source(sdf)
    directory = tmp_path / "docking"
    directory.mkdir()
    coords = ethanol.GetConformer().GetPositions()[:3] + [20, 15, 10]
    pose = directory / "protein__S1__ethanol_12345678.pdbqt"
    pose.write_text('MODEL 2\nREMARK VINA RESULT: -5 0 0\n' +
        '\n'.join(atom(i + 1, ethanol.GetAtomWithIdx(i).GetSymbol(), xyz)
                  for i, xyz in enumerate(coords)) + '\nENDMDL\n')
    config = {"system_type": "soluble", "receptor_reviewed": True, "receptor_id": "protein",
              "engine": "vina", "receptor_pdb": str(receptor),
              "docking_receptor_pdbqt": str(docked), "ligands": {
                  "ethanol": {"sdf": str(sdf), "net_charge": 0,
                              "atom_map": {"1": 0, "2": 1, "3": 2}}}}
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(config))
    return path, directory, config


def test_prepare_isolated_snapshot_preserves_scored_output(tmp_path, monkeypatch):
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    config, docking, _ = protocol_inputs(tmp_path)
    jobdir = workflow.prepare(config, docking)
    job = json.loads((jobdir / "job.json").read_text())
    assert jobdir.is_relative_to(root / "workspace")
    assert job["systems"][0]["pose"]["model"] == 2
    assert job["systems"][0]["pose"]["binding_site"] == "1"
    assert job["systems"][0]["pose"]["artifact"].startswith(str(docking))
    assert not (jobdir / "01_ethanol/state.json").exists()
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(jobdir / "01_ethanol/posed_ligand.mol"), removeHs=False)
    assert mol.GetConformer().GetPositions()[:3].mean(axis=0)[0] > 19
    for name, checksum in job["immutable_inputs"].items():
        assert digest(jobdir / name) == checksum
    # Execution sees the snapshot, even if the original chemical file disappears.
    (tmp_path / "ligand.sdf").unlink()
    assert (jobdir / job["systems"][0]["chemical_input"]["sdf"]).exists()


@pytest.mark.parametrize("change, message", [
    ({"force_field": "charmm36m"}, "Automatic GAFF2"),
    ({"water_model": "opc"}, "combination"),
    ({"temperature_k": float("nan")}, "Invalid temperature"),
    ({"salt_molar": -1}, "Invalid salt"),
    ({"receptor_reviewed": False}, "receptor_reviewed"),
    ({"dt_ps": 0.004}, "2 fs"),
    ({"padding_nm": 0.5}, "cutoff"),
])
def test_invalid_scientific_settings_fail_before_build(tmp_path, change, message):
    path, _, config = protocol_inputs(tmp_path)
    config.update(change)
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match=message):
        load_protocol(path)


def test_changed_input_stops_before_any_command(tmp_path, monkeypatch):
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    config, docking, _ = protocol_inputs(tmp_path)
    jobdir = workflow.prepare(config, docking)
    (jobdir / "01_ethanol/nvt.mdp").write_text("changed settings")
    monkeypatch.setattr(workflow, "check_dependencies", lambda *args, **kwargs: True)
    with pytest.raises(ValueError, match="Input changed"):
        workflow.run(jobdir, "npt", False)
    with pytest.raises(ValueError, match="equilibration-reviewed"):
        workflow.run(jobdir, "production", False)


def test_failed_md_system_does_not_cancel_other_prepared_systems(tmp_path, monkeypatch):
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    path, directory, config = protocol_inputs(tmp_path)
    original = next(directory.glob("*.pdbqt"))
    second = directory / "protein__S1__second_12345678.pdbqt"
    second.write_text(original.read_text())
    config["ligands"]["second"] = dict(config["ligands"]["ethanol"])
    path.write_text(json.dumps(config))
    jobdir = workflow.prepare(path, directory)
    monkeypatch.setattr(workflow, "check_dependencies", lambda *a, **kw: True)
    attempted = []
    def build(runner, job, entry):
        attempted.append(entry["directory"])
        if entry["directory"].endswith("ethanol"):
            raise RuntimeError("parameterization failed")
    monkeypatch.setattr(workflow, "build_system", build)
    with pytest.raises(RuntimeError, match="1 failed systems and 1 completed"):
        workflow.run(jobdir, "build", False)
    assert attempted == ["01_ethanol", "02_second"]
    summary = json.loads((jobdir / "run_summary.json").read_text())
    assert summary["completed"] == ["02_second"]
    assert summary["failed"] == [{"system": "01_ethanol", "error": "parameterization failed"}]


@pytest.mark.parametrize("fault", ["map", "missing-entry", "missing-file", "charge"])
def test_bad_chemical_input_is_reported_without_cancelling_other_ligands(tmp_path, monkeypatch, fault):
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    path, docking, config = protocol_inputs(tmp_path)
    original = next(docking.glob("*.pdbqt"))
    (docking / "protein__S1__bad_12345678.pdbqt").write_text(original.read_text())
    entry = dict(config["ligands"]["ethanol"])
    if fault == "map":
        entry["atom_map"] = {}
    elif fault == "missing-file":
        entry["sdf"] = "missing.sdf"
    elif fault == "charge":
        entry["net_charge"] = "unknown"
    if fault != "missing-entry":
        config["ligands"]["bad"] = entry
    path.write_text(json.dumps(config))
    handoff = tmp_path / "handoff.json"
    jobdir = workflow.start(path, docking, through="prepare", result_file=handoff)
    job = json.loads((jobdir / "job.json").read_text())
    assert [entry["pose"]["ligand"] for entry in job["systems"]] == ["ethanol"]
    failure, = job["preparation_failures"]
    assert failure["ligand"] == "bad" and failure["error"]
    assert failure["pose"]["binding_site"] == "1"
    assert json.loads(handoff.read_text())["partial"] is True
    import csv
    with (jobdir / "preparation_failures.csv").open() as handle:
        row, = list(csv.DictReader(handle))
    assert row["ligand"] == "bad" and row["receptor"] == "protein"
    attempted = []
    monkeypatch.setattr(workflow, "check_dependencies", lambda *a, **kw: True)
    monkeypatch.setattr(workflow, "build_system", lambda runner, job, entry: attempted.append(entry["directory"]))
    workflow.run(jobdir, "build", False)
    assert attempted == ["01_ethanol"]
    assert json.loads((jobdir / "run_summary.json").read_text())["preparation_failures"] == [failure]


def test_all_bad_chemical_inputs_leave_failure_report(tmp_path, monkeypatch):
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    path, docking, config = protocol_inputs(tmp_path)
    config["ligands"]["ethanol"]["atom_map"] = {}
    path.write_text(json.dumps(config))
    output = root / "failed-job"
    with pytest.raises(ValueError, match="No selected ligands passed"):
        workflow.prepare(path, docking, output=output)
    assert not (output / "job.json").exists()
    failure, = json.loads((output / "preparation_failures.json").read_text())
    assert failure["ligand"] == "ethanol" and "atom_map" in failure["error"]


@pytest.mark.skipif(os.environ.get("ULTIDOCK_MD_NATIVE_TEST") != "1",
                    reason="Opt-in test requires actual GROMACS and AmberTools/ACPYPE")
def test_real_ambertools_gromacs_through_short_production(tmp_path, monkeypatch):
    """Synthetic ethanol/villin handoff: software exercise, not a binding model."""
    assert all(shutil.which(t) for t in ("gmx", "acpype", "obabel", "antechamber", "sqm"))
    path, directory, config = protocol_inputs(tmp_path)
    protein = Path(__file__).with_name("data") / "villin.pdb"
    shutil.copy2(protein, tmp_path / "protein.pdb")
    # This is a synthetic docking container in the protein coordinate frame.
    # Heavy-atom identities are copied exactly; partial charges are irrelevant here.
    lines = [line[:76] + " C" for line in protein.read_text().splitlines() if line.startswith("ATOM")]
    (tmp_path / "protein.pdbqt").write_text("\n".join(lines) + "\n")
    xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    for line in lines])
    from rdkit import Chem
    ethanol = Chem.SDMolSupplier(str(tmp_path / "ligand.sdf"), removeHs=False)[0]
    ligand = ethanol.GetConformer().GetPositions()[:3]
    ligand += [xyz[:, 0].max() + 4 - ligand[:, 0].min(), xyz[:, 1].mean(), xyz[:, 2].mean()]
    pose = directory / "protein__S1__ethanol_12345678.pdbqt"
    pose.write_text('MODEL 1\nREMARK VINA RESULT: -1 0 0\n' +
        '\n'.join(atom(i + 1, ethanol.GetAtomWithIdx(i).GetSymbol(), pos)
                  for i, pos in enumerate(ligand)) + '\nENDMDL\n')
    config.update(production_ns=0.002, nvt_ps=2, npt_ps=2, threads=1)
    path.write_text(json.dumps(config))
    root = tmp_path / "md-simulation"
    root.mkdir()
    monkeypatch.setattr(workflow, "ROOT", root)
    jobdir = workflow.start(path, directory, through="npt")
    workflow.run(jobdir, "production", True)
    state = json.loads((jobdir / "01_ethanol/state.json").read_text())
    assert all(state[s]["status"] == "complete" for s in ("build", "em", "nvt", "npt", "production"))
    assert (jobdir / "01_ethanol/production.cpt").exists()
