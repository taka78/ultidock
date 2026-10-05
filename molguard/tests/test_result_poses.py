from pathlib import Path
import queue
import runpy
import sys
from types import ModuleType

import pytest

def atom(x):
    return f"ATOM      1  C   LIG A   1    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  0.00  0.00     0.000  C"


@pytest.fixture
def worker(monkeypatch):
    from docking import config
    monkeypatch.setitem(sys.modules, "config", config)
    grids = ModuleType("make_grids")
    for name in ("HotspotGPFGenerator", "autogenerate_centers_tsv", "ensure_grids_multi_centers",
                 "ensure_grids", "ensure_whole_protein_maps"):
        setattr(grids, name, lambda *args, **kwargs: None)
    monkeypatch.setitem(sys.modules, "make_grids", grids)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "docking"))
    monkeypatch.setattr("shutil.which", lambda _: None)
    module = runpy.run_path(str(Path(__file__).resolve().parents[2] / "docking/dock_v02.py"))
    cls = module["ProcessFileThread"]
    return cls.__new__(cls), module


def gpu_files(tmp_path):
    xml = tmp_path / "receptor__S3__ligand_1234.xml"
    xml.write_text('<autodock_gpu><ligand>/inputs/ligand.pdbqt</ligand><runs>'
                   '<run id="19"><free_NRG_binding>-7.15</free_NRG_binding></run>'
                   '<run id="32"><free_NRG_binding>-7.13</free_NRG_binding></run>'
                   '</runs></autodock_gpu>')
    xml.with_name(xml.stem + "-best.pdbqt").write_text(atom(2.0) + "\n")
    xml.with_suffix(".dlg").write_text("DOCKED: MODEL 19\nDOCKED: " + atom(1.0)
        + "\nDOCKED: ENDMDL\nDOCKED: MODEL 32\nDOCKED: " + atom(2.0) + "\nDOCKED: ENDMDL\n")
    return xml


def test_gpu_score_matches_saved_coordinates_not_first_xml_run(tmp_path, worker):
    parser, _ = worker
    xml = gpu_files(tmp_path)
    rows = parser.parse_adgpu_xml(xml, "receptor", "/inputs/ligand.pdbqt")
    assert len(rows) == 2
    best = [r for r in rows if r[7]]
    assert len(best) == 1
    assert best[0][0] == "ligand-Model32-receptor"
    assert best[0][1] == -7.13
    assert best[0][4] == str(xml.with_name(xml.stem + "-best.pdbqt"))
    assert best[0][5:7] == ("3", "/inputs/ligand.pdbqt")
    assert rows[0][4] == str(xml.with_suffix(".dlg"))


@pytest.mark.parametrize("failure", ["missing", "mismatch", "ambiguous"])
def test_gpu_never_substitutes_input_for_missing_or_unmatched_pose(tmp_path, worker, failure):
    xml = gpu_files(tmp_path)
    best = xml.with_name(xml.stem + "-best.pdbqt")
    if failure == "missing":
        best.unlink()
    elif failure == "mismatch":
        best.write_text(atom(99.0))
    else:
        dlg = xml.with_suffix(".dlg")
        dlg.write_text(dlg.read_text() + "DOCKED: MODEL 40\nDOCKED: " + atom(2.0) + "\n")
    parser, _ = worker
    assert parser.parse_adgpu_xml(xml, "receptor", "/inputs/ligand.pdbqt") == []


def test_vina_best_pose_extract_preserves_correct_model_and_coordinates(tmp_path, worker):
    output = tmp_path / "receptor__S2__ligand_1234.pdbqt"
    first = "MODEL 1\nREMARK VINA RESULT: -7.2 0.0 0.0\n" + atom(1.0) + "\nENDMDL\n"
    second = "MODEL 2\nREMARK VINA RESULT: -8.1 1.0 2.0\n" + atom(2.0) + "\nENDMDL\n"
    output.write_text(first + second)
    parser, _ = worker
    rows = parser.parse_vina_output_file(output, "receptor", "/inputs/ligand.pdbqt")
    assert len(rows) == 2
    best = next(r for r in rows if r[7])
    assert best[0] == "ligand-Model2-receptor"
    assert best[1:4] == (-8.1, 1.0, 2.0)
    assert best[4].endswith("-best.pdbqt")
    assert Path(best[4]).read_text() == second
    assert output.read_text() == first + second


def test_database_worker_preserves_extended_records(tmp_path, monkeypatch, worker):
    from docking import db_manager
    monkeypatch.setattr(db_manager, "RESULTS_DIR", str(tmp_path))
    database = db_manager.DockingDatabaseManager()
    parser, module = worker
    rows = parser.parse_adgpu_xml(gpu_files(tmp_path), "receptor", "/inputs/ligand.pdbqt")
    processor = module["DockingProcessor"].__new__(module["DockingProcessor"])
    processor.db_manager = database
    processor.db_queue = queue.Queue()
    processor.db_queue.put(rows)
    processor.db_queue.put("SHUTDOWN")
    processor._db_worker()
    stored = database.connection.execute(
        "SELECT docking_file,ligand_file,is_best_pose FROM docking_results ORDER BY id").fetchall()
    assert [tuple(r) for r in stored] == [(r[4], r[6], r[7]) for r in rows]
    # Benchmark ligand labels still use the original input, across all runs.
    from benchmarks.dude_docking_benchmark import load_docking_rows
    by_ligand, by_site = load_docking_rows(Path(database.db_path))
    assert list(by_ligand) == [Path("/inputs/ligand.pdbqt")]
    assert len(by_ligand[Path("/inputs/ligand.pdbqt")]) == 2
    assert list(by_site) == ["3"]
    database.close()
