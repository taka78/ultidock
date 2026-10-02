"""GPU requests must not silently turn into CPU docking runs."""

import contextlib
import io
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from docking import setup


class BackendSelectionTests(unittest.TestCase):
    def setUp(self):
        environment = patch.dict(os.environ)
        environment.start()
        self.addCleanup(environment.stop)
        os.environ.pop("RUSTICL_ENABLE", None)
        self.output = io.StringIO()
        self.addCleanup(self.output.close)
        redirect = contextlib.redirect_stdout(self.output)
        redirect.__enter__()
        self.addCleanup(redirect.__exit__, None, None, None)

    def test_nvidia_device_selects_cuda(self):
        with patch.object(setup.shutil, "which", return_value="/usr/bin/nvidia-smi"), \
             patch.object(setup.subprocess, "run", return_value=SimpleNamespace(
                 stdout="GPU 0: Example GPU (UUID: GPU-example)\n"
             )) as run:
            self.assertEqual(setup.detect_gpu(allow_cpu_fallback=False), "CUDA")
        run.assert_called_once()
        self.assertEqual(run.call_args.args[0], ["nvidia-smi", "-L"])

    def test_opencl_gpu_does_not_require_rocm_smi(self):
        with patch.object(setup.shutil, "which",
                          side_effect=lambda name: "/usr/bin/clinfo" if name == "clinfo" else None), \
             patch.object(setup.subprocess, "run", return_value=SimpleNamespace(
                 stdout="  Device Name        Example Radeon\n  Device Type        GPU\n"
             )):
            self.assertEqual(setup.detect_gpu(allow_cpu_fallback=False), "OPENCL")

    def test_failed_nvidia_probe_still_checks_opencl(self):
        with patch.object(setup.shutil, "which", side_effect=lambda name: f"/usr/bin/{name}"), \
             patch.object(setup.subprocess, "run", side_effect=[
                 subprocess.CalledProcessError(1, ["nvidia-smi", "-L"]),
                 SimpleNamespace(stdout="  Device Type       GPU\n"),
             ]):
            self.assertEqual(setup.detect_gpu(allow_cpu_fallback=False), "OPENCL")

    def test_empty_nvidia_and_cpu_only_opencl_do_not_count_as_gpu(self):
        with patch.object(setup.shutil, "which", side_effect=lambda name: f"/usr/bin/{name}"), \
             patch.object(setup.subprocess, "run", side_effect=[
                 SimpleNamespace(stdout="No devices were found\n"),
                 SimpleNamespace(stdout="  Device Type       CPU\n"),
             ]):
            with self.assertRaisesRegex(RuntimeError, "--mode gpu requested"):
                setup.detect_gpu(allow_cpu_fallback=False)

    def test_auto_can_fall_back_to_cpu(self):
        with patch.object(setup.shutil, "which", return_value=None):
            self.assertEqual(setup.detect_gpu(), "CPU")

    def test_gpu_setup_stops_before_build_or_config_when_no_gpu_is_detected(self):
        with patch.object(setup.shutil, "which", return_value=None), \
             patch.object(setup, "detect_and_compile_autodock_gpu") as build, \
             patch.object(setup, "_resolve_directory") as directories:
            with self.assertRaisesRegex(RuntimeError, "--mode gpu requested"):
                setup.run_setup(SimpleNamespace(mode="gpu"))
            build.assert_not_called()
            directories.assert_not_called()

    def test_explicit_backends_skip_automatic_detection(self):
        # Stop immediately after selection, before any filesystem/build work.
        for mode in ("cpu", "cuda", "opencl"):
            with self.subTest(mode=mode), \
                 patch.object(setup, "detect_gpu") as detect, \
                 patch.object(setup, "detect_opencl_gpu") as opencl, \
                 patch.object(setup, "_resolve_directory", side_effect=StopIteration):
                with self.assertRaises(StopIteration):
                    setup.run_setup(SimpleNamespace(mode=mode))
                detect.assert_not_called()
                self.assertEqual(opencl.call_count, int(mode == "opencl"))

    def test_rusticl_retry_enables_driver_only_after_gpu_appears(self):
        with patch.object(setup.shutil, "which", return_value="/usr/bin/clinfo"), \
             patch.object(setup.subprocess, "run", side_effect=[
                 SimpleNamespace(stdout="  Platform Name  rusticl\nNumber of devices 0\n"),
                 SimpleNamespace(stdout="  Device Type  GPU\n"),
             ]) as run:
            self.assertTrue(setup.detect_opencl_gpu())
            self.assertEqual(os.environ["RUSTICL_ENABLE"], "radeonsi")
            self.assertNotIn("RUSTICL_ENABLE", run.call_args_list[0].kwargs["env"])
            self.assertEqual(run.call_args_list[1].kwargs["env"]["RUSTICL_ENABLE"], "radeonsi")

    def test_platform_without_gpu_never_selects_opencl_or_keeps_retry_setting(self):
        with patch.object(setup.shutil, "which",
                          side_effect=lambda name: "/usr/bin/clinfo" if name == "clinfo" else None), \
             patch.object(setup.subprocess, "run", return_value=SimpleNamespace(
                 stdout="  Platform Name  rusticl\nNumber of devices 0\n"
             )) as run:
            with self.assertRaisesRegex(RuntimeError, "--mode gpu requested"):
                setup.detect_gpu(allow_cpu_fallback=False)
            self.assertEqual(run.call_count, 2)
            self.assertNotIn("RUSTICL_ENABLE", os.environ)

    def test_failed_retry_does_not_change_environment(self):
        with patch.object(setup.shutil, "which", return_value="/usr/bin/clinfo"), \
             patch.object(setup.subprocess, "run", side_effect=[
                 SimpleNamespace(stdout="  Platform Name  rusticl\nNumber of devices 0\n"),
                 subprocess.CalledProcessError(1, ["clinfo", "--human"]),
             ]):
            self.assertFalse(setup.detect_opencl_gpu())
            self.assertNotIn("RUSTICL_ENABLE", os.environ)

    def test_explicit_rusticl_setting_is_never_overridden(self):
        for value in ("iris", ""):
            with self.subTest(value=value), patch.dict(os.environ, RUSTICL_ENABLE=value), \
                 patch.object(setup.shutil, "which", return_value="/usr/bin/clinfo"), \
                 patch.object(setup.subprocess, "run", return_value=SimpleNamespace(
                     stdout="  Platform Name  rusticl\nNumber of devices 0\n"
                 )) as run:
                self.assertFalse(setup.detect_opencl_gpu())
                run.assert_called_once()
                self.assertEqual(os.environ["RUSTICL_ENABLE"], value)

    def test_setup_saves_runtime_setting_and_config_restores_it(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / "docking").mkdir()
            args = setup.build_parser().parse_args([
                "--mode", "gpu", "--skip-wget", "--skip-profile",
            ])
            with patch.object(setup, "ROOT_DIR", directory), \
                 patch.object(setup, "_resolve_directory", return_value=directory), \
                 patch.object(setup, "detect_and_compile_autodock_gpu") as build, \
                 patch.object(setup, "check_and_fix_receptors"), \
                 patch.object(setup.shutil, "which",
                              side_effect=lambda name: "/usr/bin/clinfo" if name == "clinfo" else None), \
                 patch.object(setup.subprocess, "run", side_effect=[
                     SimpleNamespace(stdout="  Platform Name  rusticl\nNumber of devices 0\n"),
                     SimpleNamespace(stdout="  Device Type  GPU\n"),
                 ]):
                # The build and subsequent children see the successful runtime setting.
                build.side_effect = lambda *args: self.assertEqual(
                    os.environ.get("RUSTICL_ENABLE"), "radeonsi"
                )
                values = setup.run_setup(args)
            self.assertEqual(values["GPU_TYPE"], "OPENCL")
            self.assertEqual(values["RUSTICL_ENABLE"], "radeonsi")
            config = Path(directory) / "docking/config.py"
            del os.environ["RUSTICL_ENABLE"]
            runpy.run_path(str(config))
            self.assertEqual(os.environ["RUSTICL_ENABLE"], "radeonsi")
            os.environ["RUSTICL_ENABLE"] = "iris"
            runpy.run_path(str(config))
            self.assertEqual(os.environ["RUSTICL_ENABLE"], "iris")


if __name__ == "__main__":
    unittest.main()
