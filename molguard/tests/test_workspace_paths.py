"""Workflows use portable directories and preserve existing run data."""

import pytest

from ultidock import paths


def _checkout(root):
    (root / "docking").mkdir(parents=True)
    (root / "pyproject.toml").touch()
    (root / "docking/run.py").touch()
    return root


@pytest.mark.parametrize("subdirectory", [".", "docking", "examples/nested"])
def test_regular_install_discovers_working_checkout(monkeypatch, tmp_path, subdirectory):
    checkout = _checkout(tmp_path / "checkout")
    cwd = checkout / subdirectory
    cwd.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(paths, "__file__", str(tmp_path / "site-packages/ultidock/paths.py"))
    monkeypatch.delenv("ULTIDOCK_HOME", raising=False)
    monkeypatch.chdir(cwd)

    assert paths.source_root() == checkout
    assert paths.application_root() == checkout
    assert not (checkout / ".ultidock-runtime.json").exists()


@pytest.mark.parametrize(
    "package_path", ["ultidock/paths.py", ".venv/lib/site-packages/ultidock/paths.py"]
)
def test_package_ancestors_locate_checkout_from_elsewhere(monkeypatch, tmp_path, package_path):
    checkout = _checkout(tmp_path / "checkout")
    other = _checkout(tmp_path / "other")
    monkeypatch.setattr(paths, "__file__", str(checkout / package_path))
    monkeypatch.delenv("ULTIDOCK_HOME", raising=False)
    monkeypatch.chdir(other)

    assert paths.workspace_root() == checkout


@pytest.mark.parametrize("marker", ["pyproject.toml", "docking/run.py"])
def test_incomplete_checkout_uses_current_directory(monkeypatch, tmp_path, marker):
    partial = tmp_path / marker
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.touch()
    cwd = tmp_path / "work"
    cwd.mkdir()
    monkeypatch.setattr(paths, "__file__", str(tmp_path / "site-packages/ultidock/paths.py"))
    monkeypatch.delenv("ULTIDOCK_HOME", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "data"))
    monkeypatch.chdir(cwd)

    assert paths.source_root() is None
    assert paths.workspace_root() == cwd
    assert not (tmp_path / "data").exists()


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
