"""Installed workflows must write to their managed home, preserving run data."""

from ultidock import paths


def test_installed_home_uses_xdg_data_directory(monkeypatch, tmp_path):
    monkeypatch.delenv("ULTIDOCK_HOME", raising=False)
    monkeypatch.setattr(paths, "source_root", lambda: None)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "data"))

    assert paths.workspace_root() == tmp_path / "data/ultidock"


def test_explicit_home_materializes_resources_and_preserves_results(monkeypatch, tmp_path):
    checkout = tmp_path / "checkout"
    (checkout / "docking").mkdir(parents=True)
    (checkout / "docking/run.py").write_text("# original workflow\n")
    (checkout / "docking/config.py").write_text("# unrelated source configuration\n")
    native = checkout / "docking/AUTODOCK_GPU_DIR/autogrid/ad4_shared/paramdat2h.csh"
    native.parent.mkdir(parents=True)
    native.write_text("# native parameter generator\n")
    home = tmp_path / "home"
    monkeypatch.setattr(paths, "source_root", lambda: checkout)
    monkeypatch.setenv("ULTIDOCK_HOME", str(home))

    assert paths.application_root() == home
    assert (home / native.relative_to(checkout)).is_file()
    assert not (home / "docking/config.py").exists()
    config = home / "docking/config.py"
    config.write_text("# saved configuration\n")
    result = home / "docking/RESULTS_DIR/poses.csv"
    result.parent.mkdir()
    result.write_text("saved poses\n")
    script = home / "docking/run.py"
    original_mtime = script.stat().st_mtime_ns

    paths.application_root()
    assert script.stat().st_mtime_ns == original_mtime
    (checkout / "docking/run.py").write_text("# updated workflow\n")
    paths.application_root()

    assert script.read_text() == "# updated workflow\n"
    assert config.read_text() == "# saved configuration\n"
    assert result.read_text() == "saved poses\n"


def test_path_programs_are_linked_without_replacing_local_tools(monkeypatch, tmp_path):
    programs = tmp_path / "programs"
    programs.mkdir()
    vina = programs / "vina"
    vina.write_text("#!/bin/sh\nexit 0\n")
    vina.chmod(0o755)
    home = tmp_path / "home"
    monkeypatch.setattr(paths.shutil, "which", lambda name: str(vina) if name == "vina" else None)

    paths._link_native_programs(home, None)
    linked = home / "docking/VINA_DIR/bin/vina"
    assert linked.resolve() == vina
    linked.unlink()
    linked.write_text("local installation\n")
    paths._link_native_programs(home, None)
    assert linked.read_text() == "local installation\n"
