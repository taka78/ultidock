import argparse
import importlib
import json
import os
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

if __package__:
    from .setup import build_parser as build_setup_parser, run_setup
    from .pose_provenance import record_docking_run
else:
    from setup import build_parser as build_setup_parser, run_setup
    from pose_provenance import record_docking_run

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
    parser.add_argument(
        "--keep-artifacts",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Keep original ligand archives (default for docking runs)",
    )
    add_md_arguments(parser)
    return parser


def add_md_arguments(parser):
    """Shared continuation flags for the pipeline and bundled docking examples."""
    md = parser.add_argument_group("Docking to molecular dynamics")
    md.add_argument("--md-config", "--md", type=Path, metavar="PROTOCOL",
                    help="Continue docking into GROMACS using this MD protocol JSON")
    md.add_argument("--md-through", choices=["prepare", "build", "em", "nvt", "npt"],
                    help="Last MD stage for this new job (default: npt)")
    md.add_argument("--md-work-dir", type=Path, help="New work directory under md-simulation/")
    for tool in ("gmx", "acpype", "obabel"):
        md.add_argument(f"--md-{tool}", help=f"MD {tool} executable override")


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


def _md_command(args, command, *extra):
    """Keep MD implementation in md-simulation and propagate its exit status."""
    script = SCRIPT_DIR.parent / "md-simulation" / "workflow.py"
    argv = [sys.executable, "-u", str(script), command, "--config", str(args.md_config)]
    if args.md_work_dir:
        argv.extend(["--work-dir", str(args.md_work_dir)])
    for tool in ("gmx", "acpype", "obabel"):
        value = getattr(args, f"md_{tool}")
        if value:
            argv.extend([f"--{tool}", value])
    subprocess.run([*argv, *map(str, extra)], check=True)


def main(argv: list[str] | None = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.md_config is None and any((args.md_through, args.md_work_dir, args.md_gmx,
                                       args.md_acpype, args.md_obabel)):
        parser.error("MD options require --md-config PROTOCOL (or --md PROTOCOL)")
    md_error = None
    if args.md_config:
        args.md_config = args.md_config.expanduser().resolve()
        print("[pipeline] Checking MD protocol and dependencies before docking", flush=True)
        try:
            _md_command(args, "check", *(["--inputs-only"] if args.md_through == "prepare" else []))
        except subprocess.CalledProcessError:
            md_error = "MD preflight failed; see the diagnostics above. Docking continued."
            print(f"[WARN] {md_error}", flush=True)

    config = _ensure_config(args.skip_setup, args)

    if not args.skip_extract:
        from extract import main as extract_main

        keep_artifacts = (
            args.keep_artifacts if args.keep_artifacts is not None else not args.benchmark
        )
        extract_main(keep_artifacts=keep_artifacts)

    manifest = None
    if not args.skip_docking:
        if args.md_config and md_error is None:
            engine = "adgpu" if config.GPU_TYPE.upper() in {"CUDA", "NVIDIA", "OPENCL", "AMD"} else "vina"
            try:
                _md_command(args, "check", "--inputs-only", "--engine", engine)
            except subprocess.CalledProcessError:
                md_error = "MD protocol does not match the docking backend; docking continued."
                print(f"[WARN] {md_error}", flush=True)
        from dock_v02 import DockingProcessor

        processor = DockingProcessor()
        outputs = processor.run()
        manifest = record_docking_run(config.DOCKING_DIR, outputs, getattr(processor, "failures", []))
        print(f"Docking run manifest: {manifest}", flush=True)
        if not outputs:
            print("[WARN] No successful docking outputs; analysis and MD are skipped.", flush=True)
            if args.md_config:
                _record_md_skip(manifest.with_suffix(".md.json"), "No successful docking outputs")
            return 1

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

    if args.md_config:
        result_file = (manifest.with_suffix(".md.json") if manifest else
                       Path(config.DOCKING_DIR) / "md-handoff.json")
        if md_error:
            _record_md_skip(result_file, md_error)
            return 0
        if manifest:
            protocol = json.loads(args.md_config.read_text())
            records = json.loads(manifest.read_text())["outputs"]
            if not any(record.get("receptor_id") == protocol["receptor_id"] and
                       record.get("engine") == protocol["engine"] for record in records):
                _record_md_skip(result_file, "No successful docking outputs for the requested MD receptor/engine")
                return 0
        print("[pipeline] Successful docking results → pose selection → GROMACS", flush=True)
        _md_command(args, "run", "--docking-dir", config.DOCKING_DIR,
                    "--through", args.md_through or "npt", "--result-file", result_file,
                    *(["--receptor-dir", config.MACRO_MOL_DIR] if manifest else []),
                    *(["--pose-manifest", manifest] if manifest else []))
        print(f"MD handoff record: {result_file}", flush=True)


def _record_md_skip(path, reason):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"schema": 1, "status": "skipped", "reason": reason}, indent=2) + "\n")
    print(f"[WARN] MD skipped: {reason}. Handoff record: {path}", flush=True)

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode) from None
