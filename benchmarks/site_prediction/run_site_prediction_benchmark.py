#!/usr/bin/env python3
"""Download, normalize, run, evaluate, and report site-prediction benchmarks."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "benchmarks" / "results" / "site_prediction"
DEFAULT_RAW_ROOT = REPO_ROOT / "benchmarks" / "site_prediction" / "datasets" / "raw"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.site_prediction.prepare_site_datasets import (  # noqa: E402
    ensure_dataset_source,
    prepare_dataset,
)
from benchmarks.site_prediction.schema import SITE_PREDICTION_METHODS, normalize_method_name  # noqa: E402
from benchmarks.site_prediction.runtime_tools import (  # noqa: E402
    DEFAULT_AUTOGRID4,
    command_is_available,
    first_command_token,
    is_executable_path,
    resolve_autogrid4_bin,
)


def _split_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _default_output_dir() -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    return DEFAULT_OUTPUT_ROOT / f"site_prediction_{stamp}"


def _run(command: list[str], *, dry_run: bool) -> None:
    printable = shlex.join(command)
    if dry_run:
        print(f"[dry-run] {printable}")
        return
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def _resolve_autogrid4_bin(args: argparse.Namespace) -> bool:
    configured = Path(str(args.autogrid4_bin)).expanduser()
    if args.dry_run and configured.resolve() == DEFAULT_AUTOGRID4.resolve():
        if is_executable_path(DEFAULT_AUTOGRID4):
            args.autogrid4_bin = str(DEFAULT_AUTOGRID4.resolve())
        elif args.auto_build_autogrid:
            print(f"[dry-run] Would build/check bundled AutoGrid source at {DEFAULT_AUTOGRID4}")
            args.autogrid4_bin = str(DEFAULT_AUTOGRID4)
        else:
            return False
        return True

    try:
        resolved = resolve_autogrid4_bin(
            args.autogrid4_bin,
            auto_build_bundled=bool(args.auto_build_autogrid),
            numwi=int(args.numwi),
        )
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"[preflight] AutoGrid build/check failed: {exc}", file=sys.stderr)
        return False
    if resolved:
        args.autogrid4_bin = resolved
        return True
    return False


def _preflight_methods(args: argparse.Namespace, methods: list[str]) -> list[str]:
    """Validate external method executables before dataset prep starts."""
    missing_by_method: dict[str, str] = {}
    for method in methods:
        if method == "cav-emps" and not _resolve_autogrid4_bin(args):
            missing_by_method[method] = "autogrid4"
        elif method == "fpocket" and not command_is_available(args.fpocket_cmd):
            missing_by_method[method] = first_command_token(args.fpocket_cmd) or "fpocket"
        elif method == "p2rank" and not command_is_available(args.p2rank_cmd):
            missing_by_method[method] = first_command_token(args.p2rank_cmd) or "prank"

    if not missing_by_method:
        return methods

    lines = ["Missing benchmark executable(s):"]
    for method, executable in missing_by_method.items():
        lines.append(f"  - {method}: {executable}")
    lines.extend(
        [
            "Install the missing tools, pass an explicit command/path "
            "(--autogrid4-bin, --fpocket-cmd, --p2rank-cmd), or run only installed methods.",
            "For cav-emps, the bundled AutoGrid source is built automatically unless "
            "--no-auto-build-autogrid is set.",
            "Use --allow-missing-tools only when you intentionally want a partial benchmark.",
        ]
    )

    if not args.allow_missing_tools:
        raise SystemExit("\n".join(lines))

    for method, executable in missing_by_method.items():
        print(f"[skip] {method}: missing executable {executable}")
    runnable = [method for method in methods if method not in missing_by_method]
    if not runnable:
        raise SystemExit("\n".join(lines + ["No selected benchmark methods are runnable."]))
    return runnable


def _combine_prediction_files(paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    n_rows = 0
    with output_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = None
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(f"prediction file was not produced: {path}")
            with path.open("r", encoding="utf-8", newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if not reader.fieldnames:
                    raise ValueError(f"prediction file has no header: {path}")
                if writer is None:
                    writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
                if not wrote_header:
                    writer.writeheader()
                    wrote_header = True
                for row in reader:
                    writer.writerow(row)
                    n_rows += 1
    if not wrote_header or n_rows == 0:
        raise RuntimeError(
            f"No prediction rows were produced; refusing to evaluate empty predictions at {output_path}"
        )


def build_predict_command(
    *,
    dataset: str,
    method: str,
    normalized_root: Path,
    prediction_root: Path,
    args: argparse.Namespace,
) -> list[str]:
    command = [
        sys.executable,
        str(REPO_ROOT / "benchmarks" / "site_prediction" / "predict_sites.py"),
        "--dataset",
        dataset,
        "--method",
        method,
        "--normalized-root",
        str(normalized_root),
        "--output-dir",
        str(prediction_root),
        "--jobs",
        str(args.jobs),
        "--autosites",
        str(args.autosites),
        "--seed",
        str(args.seed),
        "--autogrid4-bin",
        str(args.autogrid4_bin),
        "--numwi",
        str(args.numwi),
        "--hotspot-box-size",
        str(args.hotspot_box_size),
        "--grid-spacing",
        str(args.grid_spacing),
        "--blind-cap",
        str(args.blind_cap),
    ]
    if args.targets:
        command.extend(["--targets", args.targets])
    if args.max_targets is not None:
        command.extend(["--max-targets", str(args.max_targets)])
    if args.force:
        command.append("--force")
    if args.keep_artifacts:
        command.append("--keep-artifacts")
    if args.work_root is not None:
        command.extend(["--work-root", str(Path(args.work_root).resolve())])
    if not args.auto_build_autogrid:
        command.append("--no-auto-build-autogrid")
    if args.r_min is not None:
        command.extend(["--r-min", str(args.r_min)])
    if args.receptor_prepare_command:
        command.extend(["--receptor-prepare-command", args.receptor_prepare_command])
    if args.fpocket_cmd:
        command.extend(["--fpocket-cmd", args.fpocket_cmd])
    if args.p2rank_cmd:
        command.extend(["--p2rank-cmd", args.p2rank_cmd])
    return command


def run_site_prediction_benchmark(args: argparse.Namespace) -> Path:
    output_dir = Path(args.output_dir).resolve() if args.output_dir else _default_output_dir().resolve()
    raw_root = Path(args.raw_root).resolve()
    normalized_root = Path(args.normalized_root).resolve() if args.normalized_root else output_dir / "normalized"
    manifest_root = output_dir / "manifests"
    prediction_root = output_dir / "predictions"
    evaluation_root = output_dir / "evaluation"
    report_root = output_dir / "reports"

    datasets = _split_csv(args.datasets)
    methods = [normalize_method_name(method) for method in _split_csv(args.methods)]
    unsupported = [method for method in methods if method not in SITE_PREDICTION_METHODS]
    if unsupported:
        raise SystemExit(f"Unsupported method(s): {', '.join(unsupported)}")
    methods = _preflight_methods(args, methods)

    print(f"Output dir: {output_dir}")
    print(f"Datasets:   {', '.join(datasets)}")
    print(f"Methods:    {', '.join(methods)}")

    source_root = ensure_dataset_source(
        source_root=Path(args.source_root).resolve() if args.source_root else None,
        raw_root=raw_root,
        download=bool(args.download),
        dry_run=bool(args.dry_run),
    )
    if args.dry_run and not source_root.exists():
        print(f"[dry-run] Would download rdk/p2rank-datasets to {source_root}")

    for dataset in datasets:
        print(f"\n[{dataset}] preparing dataset")
        if source_root.exists():
            prepare_dataset(
                dataset=dataset,
                source_root=source_root,
                normalized_root=normalized_root,
                manifest_root=manifest_root,
                max_targets=args.max_targets,
                dry_run=bool(args.dry_run),
            )
        elif not args.dry_run:
            raise FileNotFoundError(f"Dataset source not found: {source_root}")

        prediction_files: list[Path] = []
        for method in methods:
            print(f"\n[{dataset}] running {method}")
            command = build_predict_command(
                dataset=dataset,
                method=method,
                normalized_root=normalized_root,
                prediction_root=prediction_root / dataset,
                args=args,
            )
            _run(command, dry_run=bool(args.dry_run))
            prediction_files.append(prediction_root / dataset / method / "predictions.tsv")

        combined_predictions = evaluation_root / dataset / "predictions.tsv"
        if args.dry_run:
            print(f"[dry-run] Would combine predictions into {combined_predictions}")
        else:
            _combine_prediction_files(prediction_files, combined_predictions)

        print(f"\n[{dataset}] evaluating")
        _run(
            [
                sys.executable,
                str(REPO_ROOT / "benchmarks" / "site_prediction" / "evaluate_predictions.py"),
                "--dataset",
                dataset,
                "--normalized-root",
                str(normalized_root),
                "--predictions-tsv",
                str(combined_predictions),
                "--output-dir",
                str(evaluation_root / dataset),
                "--threshold",
                str(args.threshold),
            ],
            dry_run=bool(args.dry_run),
        )

        print(f"\n[{dataset}] reporting")
        _run(
            [
                sys.executable,
                str(REPO_ROOT / "benchmarks" / "site_prediction" / "report.py"),
                "--evaluation-dir",
                str(evaluation_root / dataset),
                "--output-dir",
                str(report_root / dataset),
            ],
            dry_run=bool(args.dry_run),
        )

    return output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the automatic COACH420/HOLO4K site-prediction benchmark.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--datasets", default="coach420,holo4k", help="Comma-separated datasets.")
    parser.add_argument("--methods", default="cav-emps,fpocket,p2rank", help="Comma-separated methods.")
    parser.add_argument("--source-root", type=Path, help="Existing rdk/p2rank-datasets checkout.")
    parser.add_argument("--raw-root", default=str(DEFAULT_RAW_ROOT), type=Path)
    parser.add_argument("--normalized-root", type=Path, help="Override normalized dataset output root.")
    parser.add_argument("--output-dir", type=Path, help="Benchmark output directory.")
    parser.add_argument("--download", dest="download", action="store_true", default=True)
    parser.add_argument("--no-download", dest="download", action="store_false")
    parser.add_argument("--targets", help="Comma-separated target IDs or all.")
    parser.add_argument("--max-targets", type=int, help="Limit targets per dataset for smoke tests.")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--autosites", type=int, default=6)
    parser.add_argument("--threshold", type=float, default=4.0, help="DCA hit threshold in Angstrom.")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--keep-artifacts",
        action="store_true",
        help=(
            "Keep per-target CaV-EMPS AutoGrid maps in the prediction result tree. "
            "By default, heavy map artifacts are written to temporary workdirs and removed."
        ),
    )
    parser.add_argument(
        "--full-autogrid-maps",
        action="store_true",
        help=(
            "Deprecated compatibility flag. CaV-EMPS benchmarks always generate "
            "full AD4 ligand maps because AutoGrid dsolvmap depends on the "
            "requested ligand type set."
        ),
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        help=(
            "Temporary work root passed to method prediction jobs when --keep-artifacts is not set. "
            "Defaults inside each method prediction output directory."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--autogrid4-bin",
        default=str(REPO_ROOT / "docking" / "AUTODOCK_GPU_DIR" / "autogrid" / "autogrid4"),
    )
    parser.add_argument(
        "--no-auto-build-autogrid",
        dest="auto_build_autogrid",
        action="store_false",
        help="Do not compile the bundled AutoGrid source automatically for cav-emps.",
    )
    parser.set_defaults(auto_build_autogrid=True)
    parser.add_argument("--numwi", type=int, default=128, help="NUMWI passed to Ultidock's build script.")
    parser.add_argument("--grid-spacing", type=float, default=0.375)
    parser.add_argument("--blind-cap", type=float, default=150.0)
    parser.add_argument("--r-min", type=float)
    parser.add_argument("--hotspot-box-size", type=float, default=35.0)
    parser.add_argument("--receptor-prepare-command")
    parser.add_argument("--fpocket-cmd", default="fpocket")
    parser.add_argument("--p2rank-cmd", default="prank predict")
    parser.add_argument(
        "--allow-missing-tools",
        action="store_true",
        help="Skip selected methods whose executable is missing instead of failing before setup.",
    )
    return parser


def main() -> None:
    output_dir = run_site_prediction_benchmark(build_parser().parse_args())
    print(f"\nBenchmark output: {output_dir}")


if __name__ == "__main__":
    main()
