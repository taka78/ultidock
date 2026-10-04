"""The source files and fixtures required by the workflow scripts."""

from pathlib import Path


def runtime_files(root: Path) -> list[Path]:
    """Select resources independently of Git, also when building an sdist."""
    excluded_dirs = {
        "__pycache__",
        "workspace",
        "results",
        "datasets",
        "backups",
        "provenance",
        "expected-results",
        "tests",
        "bin",
        "autom4te.cache",
        "DOCKING_DIR",
        "RESULTS_DIR",
        "ANALYSIS_DIR",
        "LIGANDS_DIR",
        "MACRO_MOL_DIR",
    }
    generated_names = {
        "config.py",
        "config.log",
        "config.status",
        "autogrid4",
        "test_cuda",
        "default_parameters.h",
        "stringify.h",
    }
    generated_suffixes = {
        ".pyc",
        ".o",
        ".a",
        ".so",
        ".db",
        ".sqlite",
        ".log",
        ".map",
        ".fld",
        ".gpf",
        ".glg",
        ".npy",
        ".npz",
    }
    files = [
        root / name
        for name in ("README.md", "SETUP.md", "LICENSE", "CHANGELOG.md")
        if (root / name).is_file()
    ]
    for directory in (
        "ultidock",
        "cli",
        "molguard",
        "docking",
        "benchmarks",
        "examples",
        "scripts",
        "md-simulation",
    ):
        for path in (root / directory).rglob("*"):
            relative = path.relative_to(root)
            if not path.is_file():
                continue
            if path.is_symlink() and root.resolve() not in path.resolve().parents:
                continue
            if any(part in excluded_dirs or part.startswith(".") for part in relative.parts):
                continue
            if directory == "md-simulation" and "tools" in relative.parts:
                continue
            if path.name in generated_names or path.suffix in generated_suffixes:
                continue
            if relative.as_posix() == "docking/AUTODOCK_GPU_DIR/autogrid/Makefile":
                continue
            if directory == "ultidock" and path.suffix != ".py":
                continue
            if directory == "docking" and "AUTODOCK_GPU_DIR" not in relative.parts:
                if path.suffix not in {".py", ".sh"} and path.name != "LICENSE":
                    continue
            if directory in {"cli", "molguard", "benchmarks"} and path.suffix not in {".py", ".md"}:
                continue
            if directory == "md-simulation" and path.suffix not in {".py", ".md", ".json", ".mdp"}:
                continue
            if directory == "examples" and path.suffix not in {
                ".py",
                ".md",
                ".json",
                ".pdbqt",
                ".pdb",
                ".mol2",
            }:
                continue
            files.append(path)
    return sorted(files)
