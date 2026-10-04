"""Exercise the real docking parser methods without a user's runtime config."""

import ast
from pathlib import Path
import runpy
import queue
import sqlite3
import sys
from types import ModuleType, SimpleNamespace

import pytest

from test_handoff import vina


def parser(monkeypatch, tmp_path):
    root = Path(__file__).resolve().parents[2] / "docking"
    config = ModuleType("config")
    tree = ast.parse((root / "dock_v02.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "config":
            for alias in node.names:
                setattr(config, alias.name, str(tmp_path) if alias.name.endswith("DIR") else 1)
    config.GPU_TYPE, config.NUMWI = "CPU", "64"
    config.GRID_MODE = "known_site"
    monkeypatch.setitem(sys.modules, "config", config)
    grids = ModuleType("make_grids")
    for name in ("HotspotGPFGenerator", "autogenerate_centers_tsv", "ensure_grids_multi_centers",
                 "ensure_grids", "ensure_whole_protein_maps"):
        setattr(grids, name, None)
    monkeypatch.setitem(sys.modules, "make_grids", grids)
    monkeypatch.syspath_prepend(str(root))
    cls = runpy.run_path(str(root / "dock_v02.py"))["ProcessFileThread"]
    return cls.__new__(cls)


def test_vina_database_path_is_docked_output(monkeypatch, tmp_path):
    output = tmp_path / "protein__S2__ethanol_12345678.pdbqt"
    input_ligand = tmp_path / "ethanol.pdbqt"
    vina(output, [-8, -7])
    rows = parser(monkeypatch, tmp_path).parse_vina_output_file(output, "protein", input_ligand)
    assert len(rows) == 2 and rows[0][0] == "ethanol-Model1-protein"
    assert rows[0][4] == str(output.resolve()) and rows[0][4] != str(input_ligand)
    assert rows[0][5] == "2"


def test_gpu_database_path_is_dlg_coordinates(monkeypatch, tmp_path):
    xml = tmp_path / "protein__S3__ethanol_12345678.xml"
    xml.write_text('<autodock_gpu><ligand>ethanol.pdbqt</ligand><runs><run id="42">'
                   '<free_NRG_binding>-9</free_NRG_binding></run></runs></autodock_gpu>')
    rows = parser(monkeypatch, tmp_path).parse_adgpu_xml(xml, "protein", "ethanol.pdbqt")
    assert rows[0][4] == str(xml.with_suffix(".dlg").resolve())
    assert rows[0][0] == "ethanol-Model42-protein" and rows[0][5] == "3"


@pytest.mark.parametrize("failure", ["exit", "timeout", "empty"])
def test_real_docking_worker_records_case_failures(monkeypatch, tmp_path, failure):
    worker = parser(monkeypatch, tmp_path)
    namespace = worker.run.__globals__
    ligand = tmp_path / "ethanol.pdbqt"
    ligand.write_text("fixture")
    worker._parent = SimpleNamespace(receptor_sites=[(tmp_path / "protein.pdbqt", [
        {"center": [0, 0, 0], "npts": [41, 41, 41], "spacing": 0.375,
         "site_id": "S1", "fld_path": str(tmp_path / "maps.fld")}
    ])])
    worker.bunch = [str(ligand)]
    worker.db_queue = queue.Queue()
    worker.vina_sem = namespace["Semaphore"](1)
    worker.callback = lambda: None
    monkeypatch.setitem(namespace, "_normalize_ligand", lambda path: True)
    monkeypatch.setattr(worker, "parse_vina_output_file", lambda *a: [])
    def command(*args, **kwargs):
        if failure == "timeout":
            raise namespace["subprocess"].TimeoutExpired("vina", 1)
        return SimpleNamespace(returncode=1 if failure == "exit" else 0, stdout="", stderr="failed")
    monkeypatch.setattr(namespace["subprocess"], "run", command)
    assert worker.run() == []
    assert len(worker.failures) == 1
    assert worker.failures[0]["stage"] == "docking"
    assert worker.failures[0]["ligand_id"] == "ethanol"
    assert worker.failures[0]["receptor_id"] == "protein"


def test_pool_records_worker_failure_after_database_shutdown(monkeypatch, tmp_path):
    worker = parser(monkeypatch, tmp_path)
    namespace = worker.run.__globals__
    cls = namespace["DockingProcessor"]
    processor = cls.__new__(cls)
    processor.FILES = [str(tmp_path / "ligand.pdbqt")]
    processor.barrier = processor.event = processor.vina_sem = None
    processor.MACRO_MOL_DIR = str(tmp_path)
    processor.receptor_sites = [(tmp_path / "protein.pdbqt", [])]
    processor.failures = []
    processor.db_errors = []
    processor.db_queue = queue.Queue()
    processor.db_manager = SimpleNamespace(insert_bulk=lambda records: None, close=lambda: None)
    processor.db_thread = namespace["threading"].Thread(target=processor._db_worker)
    processor.db_thread.start()
    monkeypatch.setattr(processor, "check_memory", lambda: None)
    def fail(self):
        raise RuntimeError("worker failure")
    monkeypatch.setattr(type(worker), "run", fail)
    assert processor.process_files() == []
    assert processor.failures[0]["error"] == "worker failure"
    assert not processor.db_thread.is_alive()
    assert processor.db_queue.unfinished_tasks == 0


def test_worker_continues_other_sites_receptors_and_ligands(monkeypatch, tmp_path):
    worker = parser(monkeypatch, tmp_path)
    namespace = worker.run.__globals__
    worker._parent = SimpleNamespace(receptor_sites=[(tmp_path / f"{name}.pdbqt", [
        {"site_id": "S1"}, {"site_id": "S2"}]) for name in ("one", "two")])
    worker.bunch = [str(tmp_path / f"{name}.pdbqt") for name in ("bad", "good", "next")]
    worker.db_queue = queue.Queue()
    worker.callback = lambda: None
    monkeypatch.setitem(namespace, "_normalize_ligand", lambda path: path.stem != "bad")
    attempts = []
    def dock(ligand, receptor, site):
        case = (Path(ligand).stem, receptor.stem, site["site_id"])
        attempts.append(case)
        if case == ("good", "one", "S1"):
            raise RuntimeError("timeout")
        return "__".join(case), [("pose", -5, 0, 0, "out", "1")]
    monkeypatch.setattr(worker, "_dock_site", dock)
    assert len(worker.run()) == 7
    assert len(attempts) == 8 and attempts[-1] == ("next", "two", "S2")
    assert len(worker.failures) == 2
    assert worker.db_queue.qsize() == 7


def test_bad_receptor_grids_do_not_prevent_other_receptors(monkeypatch, tmp_path):
    worker = parser(monkeypatch, tmp_path)
    load_sites = worker.run.__globals__["load_receptor_sites"]
    for name in ("bad", "good"):
        (tmp_path / f"{name}.pdbqt").write_text("fixture")
    fld = tmp_path / "good.fld"
    fld.write_text("fixture")
    def load(path, **kwargs):
        if Path(path).stem == "bad":
            raise RuntimeError("grid preparation failed")
        return [{"fld_path": str(fld)}]
    failures = []
    pairs = load_sites(str(tmp_path), load, failures=failures)
    assert [path.stem for path, _ in pairs] == ["good"]
    assert failures[0]["receptor_id"] == "bad"


def test_database_failure_is_reported_without_hanging_queue(monkeypatch, tmp_path):
    worker = parser(monkeypatch, tmp_path)
    cls = worker.run.__globals__["DockingProcessor"]
    processor = cls.__new__(cls)
    processor.db_queue, processor.db_errors = queue.Queue(), []
    def fail(records):
        raise sqlite3.OperationalError("disk full")
    processor.db_manager = SimpleNamespace(insert_bulk=fail)
    processor.db_queue.put([("pose", -5, 0, 0, "output.pdbqt", "1")])
    processor.db_queue.put("SHUTDOWN")
    processor._db_worker()
    assert processor.db_errors == ["disk full"]
    assert processor.db_queue.unfinished_tasks == 0


def test_bulk_insert_raises_on_storage_failure(monkeypatch, tmp_path):
    from docking import db_manager
    monkeypatch.setattr(db_manager, "RESULTS_DIR", str(tmp_path))
    manager = db_manager.DockingDatabaseManager("test.db")
    manager.connection.execute("PRAGMA query_only = ON")
    try:
        with pytest.raises(sqlite3.OperationalError):
            manager.insert_bulk([("pose", -5, 0, 0, "output.pdbqt", "1")])
    finally:
        manager.close()
