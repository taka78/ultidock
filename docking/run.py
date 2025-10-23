import argparse
import importlib
import os
import sys
from contextlib import contextmanager
from pathlib import Path

from setup import build_parser as build_setup_parser, run_setup

SCRIPT_DIR = Path(__file__).resolve().parent


def build_parser() -> argparse.ArgumentParser:
    parser = build_setup_parser()
    parser.description = "Ultidock pipeline runner"
    parser.add_argument(
        "--skip-setup",
        action="store_true",
        help="Skip setup step (requires existing config.py)",
    )
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip ligand extraction (vina_split) stage",
    )
    parser.add_argument(
        "--skip-docking",
        action="store_true",
        help="Skip docking stage",
    )
    parser.add_argument(
        "--skip-analysis",
        action="store_true",
        help="Skip results analysis stage",
    )
    return parser


@contextmanager
def _temporary_argv(argv):
    original = sys.argv[:]
    sys.argv = argv
    try:
        yield
    finally:
        sys.argv = original


def _ensure_config(skip_setup: bool, args: argparse.Namespace):
    if not skip_setup:
        run_setup(args)
    config_path = SCRIPT_DIR / "config.py"
    if not config_path.exists():
        raise FileNotFoundError(
            f"config.py not found at {config_path}. Run setup or remove --skip-setup."
        )
    if "config" in sys.modules:
        return importlib.reload(sys.modules["config"])
    return importlib.import_module("config")


def main(argv: list[str] | None = None):
    parser = build_parser()
    args = parser.parse_args(argv)

    config = _ensure_config(args.skip_setup, args)

    if not args.skip_extract:
        from extract import main as extract_main

        extract_main()

    if not args.skip_docking:
        from dock_v02 import DockingProcessor

        DockingProcessor().run()

    if not args.skip_analysis:
        try:
            analysis = importlib.import_module("analyse_docking_results")
        except ModuleNotFoundError as exc:
            missing = exc.name or "unknown dependency"
            print(
                f"[WARN] Optional analysis dependencies missing ({missing}); "
                "skipping analysis stage."
            )
        else:
            with _temporary_argv([str(analysis.__file__)]):
                analysis.main()

if __name__ == "__main__":
    main()
