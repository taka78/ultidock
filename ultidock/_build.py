"""Bundle legacy workflow resources without local results or native binaries."""

from pathlib import Path
import runpy
from zipfile import ZIP_DEFLATED, ZipFile

from setuptools.command.build_py import build_py

# Setuptools loads cmdclass files before the project is importable.
runtime_files = runpy.run_path(str(Path(__file__).with_name("runtime_files.py")))["runtime_files"]


class BuildPy(build_py):
    def get_source_files(self):
        root = Path(__file__).resolve().parents[1]
        return sorted(
            set(super().get_source_files())
            | {path.relative_to(root).as_posix() for path in runtime_files(root)}
        )

    def find_package_modules(self, package, package_dir):
        return [
            item
            for item in super().find_package_modules(package, package_dir)
            if not (package == "docking" and item[1] == "config")
        ]

    def run(self):
        super().run()
        if self.editable_mode:
            return
        root = Path(__file__).resolve().parents[1]
        output = Path(self.build_lib) / "ultidock" / "runtime.zip"
        output.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(output, "w", compression=ZIP_DEFLATED) as archive:
            for path in runtime_files(root):
                archive.write(path, path.relative_to(root).as_posix())
