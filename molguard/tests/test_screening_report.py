import json

from cli.report import generate_report
from docking.pose_provenance import record_docking_run


def test_report_labels_partial_screening_and_keeps_failure_details(tmp_path):
    docking = tmp_path / "docking"
    docking.mkdir()
    output = docking / "successful.pdbqt"
    output.write_text("test output")
    manifest = record_docking_run(docking, [output], [{"stage": "docking", "ligand_id": "failed",
        "receptor_id": "protein", "binding_site": "2", "error": "Vina timed out"}])
    report = generate_report(tmp_path)["report_md"].read_text()
    assert "Status: partial; successful outputs: 1; failed cases: 1" in report
    assert "Vina timed out" in report and manifest.with_suffix(".failures.csv").name in report


def test_all_failed_manifest_is_durable_and_empty(tmp_path):
    path = record_docking_run(tmp_path, [], [{"stage": "receptor-grids",
        "receptor_id": "protein", "error": "missing maps"}])
    data = json.loads(path.read_text())
    assert data["status"] == "failed" and data["outputs"] == []
    assert data["failed_cases"] == 1
    assert "missing maps" in path.with_suffix(".failures.csv").read_text()
