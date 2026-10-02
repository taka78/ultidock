"""Optional pocket tools must not block unrelated docking workflows."""

import contextlib
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from docking import pocket_boxes, setup


def executable(path: Path, body: str = "exit 0\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


class PocketToolInstallationTests(unittest.TestCase):
    def test_setup_checks_tools_without_starting_installer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable(root / "external/fpocket/bin/fpocket")
            output = io.StringIO()
            with patch.object(setup, "ROOT_DIR", directory), \
                 patch.object(setup.subprocess, "run") as run, \
                 contextlib.redirect_stdout(output):
                setup.check_pocket_tools()
            run.assert_not_called()
            self.assertIn("fpocket: ready", output.getvalue())
            self.assertIn("p2rank: not installed; will install when selected", output.getvalue())

    def test_missing_tool_installs_only_selected_method(self):
        for method in ("fpocket", "p2rank"):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                binary = Path(directory) / method
                with patch.dict(pocket_boxes.LOCAL_BINARIES, {method: binary}), \
                     patch.object(pocket_boxes.subprocess, "run",
                                  side_effect=lambda *a, **kw: executable(binary)) as run:
                    self.assertEqual(pocket_boxes.resolve_binary(method), binary)
                self.assertEqual(run.call_args.args[0][-1], method)
                run.assert_called_once()

    def test_existing_or_explicit_tool_does_not_install(self):
        with tempfile.TemporaryDirectory() as directory:
            binary = executable(Path(directory) / "tool")
            with patch.dict(pocket_boxes.LOCAL_BINARIES, {"fpocket": binary}), \
                 patch.object(pocket_boxes.subprocess, "run") as run:
                self.assertEqual(pocket_boxes.resolve_binary("fpocket"), binary)
                self.assertEqual(pocket_boxes.resolve_binary("p2rank", binary), binary)
                with self.assertRaises(FileNotFoundError):
                    pocket_boxes.resolve_binary("p2rank", binary.parent / "missing")
                run.assert_not_called()

    def test_installer_keeps_method_dependencies_separate(self):
        # Stub downloads/builds in an isolated checkout; unwanted tools fail immediately.
        installer = pocket_boxes.REPO_ROOT / "scripts/install_pocket_tools.sh"
        for method in ("fpocket", "p2rank"):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                script = root / "scripts/install_pocket_tools.sh"
                script.parent.mkdir()
                shutil.copy2(installer, script)
                stubs = root / "stubs"
                for command in ("git", "make", "curl", "tar"):
                    executable(stubs / command, "exit 99\n")
                if method == "fpocket":
                    (root / "external/fpocket").mkdir(parents=True)
                    source = root / "external/fpocket/src/fparams.c"
                    source.parent.mkdir(parents=True)
                    source.write_text(
                        "char *residue_string[M_MAX_CUSTOM_POCKET_LEN];\n"
                        "strcpy(&residue_string, pt);\n",
                        encoding="utf-8",
                    )
                    executable(stubs / "make",
                               'mkdir -p "$2/bin"\nprintf "#!/bin/bash\\nexit 0\\n" > "$2/bin/fpocket"\n'
                               'chmod +x "$2/bin/fpocket"\n')
                    expected = root / "external/fpocket/bin/fpocket"
                    unselected = root / "external/bin/prank"
                else:
                    executable(stubs / "curl")
                    executable(stubs / "tar",
                               'mkdir -p "$4/p2rank_2.5"\n'
                               'printf "#!/bin/bash\\nexit 0\\n" > "$4/p2rank_2.5/prank"\n'
                               'chmod +x "$4/p2rank_2.5/prank"\n')
                    expected = root / "external/bin/prank"
                    unselected = root / "external/fpocket"
                env = dict(os.environ, PATH=f"{stubs}:{os.environ['PATH']}")
                result = subprocess.run(["bash", str(script), method], env=env,
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(os.access(expected, os.X_OK))
                self.assertFalse(unselected.exists())
                if method == "fpocket":
                    self.assertIn("char residue_string[M_MAX_CUSTOM_POCKET_LEN];", source.read_text())
                    self.assertIn("strcpy(residue_string, pt);", source.read_text())


if __name__ == "__main__":
    unittest.main()
