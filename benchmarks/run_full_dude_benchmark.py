#!/usr/bin/env python3
"""Run the full DUD-E benchmark end-to-end without human interaction.

Usage examples:

  # Quick dev test on one target already downloaded
  python3 benchmarks/run_full_dude_benchmark.py --targets cdk2 --mode cpu

  # Download + run all 102 targets on GPU
  python3 benchmarks/run_full_dude_benchmark.py --download --mode gpu

  # Resume an interrupted run (skip targets with existing scores.csv)
  python3 benchmarks/run_full_dude_benchmark.py --download --mode gpu --continue

  # Subsample ligands for a quick feasibility check
  python3 benchmarks/run_full_dude_benchmark.py --targets cdk2,hivpr --max-ligands 50 --mode cpu

  # Dry-run: validate everything without actually docking
  python3 benchmarks/run_full_dude_benchmark.py --targets cdk2 --dry-run
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
DEFAULT_DATASET_ROOT = BENCHMARKS_DIR / "datasets"
DEFAULT_OUTPUT_ROOT = BENCHMARKS_DIR / "results" / "dude"

sys.path.insert(0, str(REPO_ROOT))

try:
    from benchmarks.download_dude import ALL_DUDE_TARGETS
except ImportError:
    from download_dude import ALL_DUDE_TARGETS

try:
    from benchmarks.summary_tools import (
        write_method_outputs,
    )
except ImportError:
    from summary_tools import (  # type: ignore
        write_method_outputs,
    )


# ── helpers ───────────────────────────────────────────────────────────────

def _check_meeko() -> tuple[str | None, str | None]:
    """Return (receptor_cmd, ligand_cmd) or (None, None) if unavailable."""
    rec_cmd = "mk_prepare_receptor.py -i {input} -p {output} --allow_bad_res"
    lig_cmd = "mk_prepare_ligand.py -i {input} -o {output}"

    # quick check: is mk_prepare_ligand.py importable / on PATH?
    if shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand"):
        return rec_cmd, lig_cmd

    # try python -m meeko.cli.mk_prepare_ligand
    try:
        subprocess.run(
            [sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-h"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        # meeko installed but scripts might not be on PATH
        rec_cmd = f"{sys.executable} -m meeko.cli.mk_prepare_receptor -i {{input}} -p {{output}} --allow_bad_res"
        lig_cmd = f"{sys.executable} -m meeko.cli.mk_prepare_ligand -i {{input}} -o {{output}}"
        return rec_cmd, lig_cmd
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    return None, None


def _check_obabel() -> tuple[str | None, str | None]:
    """Fallback: use Open Babel for conversion if Meeko is absent."""
    if shutil.which("obabel"):
        rec_cmd = "obabel {input} -O {output} -xr"
        lig_cmd = "obabel {input} -O {output}"
        return rec_cmd, lig_cmd
    return None, None


def resolve_prepare_commands(
    args: argparse.Namespace,
) -> tuple[str | None, str | None]:
    """Determine receptor/ligand prepare commands."""
    # explicit CLI overrides
    if args.receptor_prepare_command and args.ligand_prepare_command:
        return args.receptor_prepare_command, args.ligand_prepare_command

    # auto-detect
    rec, lig = _check_meeko()
    if rec:
        print("[INFO] Meeko detected — will use mk_prepare_receptor/ligand for conversion.")
        return rec, lig

    rec, lig = _check_obabel()
    if rec:
        print("[INFO] Open Babel detected — will use obabel for conversion.")
        return rec, lig

    return None, None


def discover_targets(dataset_root: Path) -> list[str]:
    """Return all target names that have required files."""
    targets: list[str] = []
    for subdir in sorted(dataset_root.iterdir()):
        if not subdir.is_dir():
            continue
        name = subdir.name.lower()
        # minimal check: at least a receptor and actives file
        has_receptor = any((subdir / f"receptor{ext}").exists() for ext in (".pdb", ".pdbqt", ".mol2"))
        has_actives = any(
            (subdir / name).exists()
            for name in ("actives_final.mol2.gz", "actives_final.pdbqt.gz", "actives.mol2.gz", "actives.pdbqt")
        )
        if has_receptor and has_actives:
            targets.append(name)
    return targets


def target_already_done(output_dir: Path, target_name: str, method_label: str) -> bool:
    """Check if a target has already been scored."""
    scores_csv = output_dir / target_name / method_label / "scores.csv"
    return scores_csv.exists()


def method_label_for_args(args: argparse.Namespace) -> str:
    """Return the result-subdirectory label for the active benchmark arm."""
    label = f"ultidock_{args.mode}"
    if getattr(args, "center_source", "auto") != "auto":
        label = f"{label}_{args.center_source}"
    return label


def write_failure_record(path: Path, payload: dict[str, object]) -> None:
    """Persist a per-target failure record for resumable benchmark runs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def move_replace(src: Path, dst: Path) -> None:
    """Move a file or directory, replacing any existing destination."""
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(str(src), str(dst))


def copy_replace(src: Path, dst: Path) -> None:
    """Copy a file, replacing any existing destination."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    shutil.copy2(str(src), str(dst))


def rewrite_metric_paths_for_method_dir(metrics: dict[str, object], target_result_dir: Path) -> None:
    """Point stored artifact paths at the method-specific output directory."""
    replacements = {
        "scores_csv": target_result_dir / "scores.csv",
        "setup_log": target_result_dir / "setup.log",
        "dock_log": target_result_dir / "dock_v02.log",
        "site_summary_json": target_result_dir / "site_summary.json",
        "site_summary_csv": target_result_dir / "site_summary.csv",
        "per_site_manifest": target_result_dir / "site_metrics.json",
        "centers_tsv": target_result_dir / "centers.tsv",
        "db_path": target_result_dir / "ultidock_results.db",
    }
    for key, path in replacements.items():
        if path.exists():
            metrics[key] = str(path.resolve())

    per_site_metrics = metrics.get("per_site_metrics")
    if not isinstance(per_site_metrics, list):
        return
    for entry in per_site_metrics:
        if not isinstance(entry, dict):
            continue
        site_id = str(entry.get("site_id") or "").strip()
        if not site_id:
            continue
        site_dir = target_result_dir / "sites" / site_id
        scores_csv = site_dir / "scores.csv"
        summary_json = site_dir / "summary.json"
        if scores_csv.exists():
            entry["scores_csv"] = str(scores_csv.resolve())
        if summary_json.exists():
            entry["summary_json"] = str(summary_json.resolve())
        db_path = target_result_dir / "ultidock_results.db"
        if db_path.exists():
            entry["db_path"] = str(db_path.resolve())


def _subsample_mol2_gz(
    original: Path, output: Path, max_records: int, seed: int = 42
) -> int:
    """Subsample a multi-record mol2.gz file to max_records.  Returns actual count."""
    import gzip
    import re

    text = gzip.open(original, "rt", encoding="utf-8").read()
    marker = "@<TRIPOS>MOLECULE"
    starts = [m.start() for m in re.finditer(re.escape(marker), text)]
    if len(starts) <= max_records:
        # just copy as-is
        shutil.copy2(original, output)
        return len(starts)

    starts.append(len(text))
    rng = random.Random(seed)
    chosen = sorted(rng.sample(range(len(starts) - 1), max_records))
    chunks = [text[starts[i] : starts[i + 1]] for i in chosen]

    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write("".join(chunks))
    return len(chunks)


def prepare_subsampled_target(
    dataset_root: Path,
    target_name: str,
    max_ligands: int,
    work_dir: Path,
    seed: int = 42,
) -> Path:
    """Create a subsampled copy of a target dataset.  Returns the subsampled target dir."""
    src = dataset_root / target_name
    dst = work_dir / target_name
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    # copy receptor and reference ligand as-is
    for fname in ("receptor.pdb", "receptor.pdbqt", "crystal_ligand.mol2",
                   "crystal_ligand.pdbqt", "benchmark.json"):
        src_file = src / fname
        if src_file.exists():
            shutil.copy2(src_file, dst / fname)

    # subsample actives
    for actives_name in ("actives_final.mol2.gz", "actives.mol2.gz"):
        src_act = src / actives_name
        if src_act.exists():
            n = _subsample_mol2_gz(src_act, dst / actives_name, max_ligands, seed)
            print(f"    Subsampled actives: {n} records")
            break

    # subsample decoys
    for decoys_name in ("decoys_final.mol2.gz", "decoys.mol2.gz"):
        src_dec = src / decoys_name
        if src_dec.exists():
            n = _subsample_mol2_gz(src_dec, dst / decoys_name, max_ligands, seed)
            print(f"    Subsampled decoys: {n} records")
            break

    return dst


# ── main pipeline ─────────────────────────────────────────────────────────

def run_benchmark_for_targets(args: argparse.Namespace) -> dict:
    """Run the full DUD-E benchmark."""
    dataset_root = Path(args.dataset_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    method_label = method_label_for_args(args)

    # ── download if requested ──
    if args.download:
        print("=" * 60)
        print("STEP 1: Downloading DUD-E datasets")
        print("=" * 60)
        dl_cmd = [
            sys.executable,
            str(BENCHMARKS_DIR / "download_dude.py"),
            "--dataset-root", str(dataset_root),
        ]
        if args.targets:
            dl_cmd.extend(["--targets", args.targets])
        subprocess.run(dl_cmd, check=True)
        print()

    # ── discover / filter targets ──
    discovered_targets = discover_targets(dataset_root)
    if args.targets and args.targets.lower() != "all":
        requested_targets = [t.strip().lower() for t in args.targets.split(",") if t.strip()]
    else:
        requested_targets = discovered_targets

    if not requested_targets:
        print("ERROR: No valid benchmark targets found.", file=sys.stderr)
        sys.exit(1)

    tracked_targets = sorted(set(discovered_targets) | set(requested_targets))
    if args.max_targets_per_run is not None and args.max_targets_per_run <= 0:
        raise ValueError("--max-targets-per-run must be a positive integer")

    print("=" * 60)
    print(f"DUD-E BENCHMARK ARM: {method_label}")
    print("=" * 60)
    for name in requested_targets:
        print(f"  • {name}")
    print()

    if args.refresh_summary:
        refreshed = write_method_outputs(
            dataset_root=dataset_root,
            output_dir=output_dir,
            method_label=method_label,
            target_names=tracked_targets,
            mode=args.mode,
            center_source=args.center_source,
            seed=args.seed,
        )
        print(f"Status:   {refreshed['paths']['status']}")
        print(f"Summary:  {refreshed['paths']['summary']}")
        print(f"Runtime:  {refreshed['paths']['runtime_summary']}")
        print(f"Failures: {refreshed['paths']['failure_summary']}")
        return refreshed["summary"]

    # ── dry-run: just validate ──
    if args.dry_run:
        # resolve prepare commands only for validation output
        rec_cmd, lig_cmd = resolve_prepare_commands(args)
        if rec_cmd is None or lig_cmd is None:
            print("WARNING: No conversion tool found (Meeko or Open Babel).", file=sys.stderr)
            print("  Receptors must already be .pdbqt and ligands must be .pdbqt", file=sys.stderr)
            print("  Install Meeko: pip install meeko", file=sys.stderr)
            print("  Or install Open Babel: apt install openbabel", file=sys.stderr)

        print("\n[DRY-RUN] Validation complete. Would process these targets:")
        for name in requested_targets:
            status = (
                "SKIP (done)"
                if getattr(args, "continue_run", False) and target_already_done(output_dir, name, method_label)
                else "PENDING"
            )
            print(f"  {name}: {status}")
        if rec_cmd:
            print(f"\n  receptor_prepare_command = {rec_cmd}")
            print(f"  ligand_prepare_command   = {lig_cmd}")
        else:
            print("\n  [WARN] No prepare commands — .pdbqt input required")
        return {"status": "dry-run", "targets": requested_targets, "method_label": method_label}

    # ── resolve prepare commands ──
    rec_cmd, lig_cmd = resolve_prepare_commands(args)
    if rec_cmd is None or lig_cmd is None:
        print("WARNING: No conversion tool found (Meeko or Open Babel).", file=sys.stderr)
        print("  Receptors must already be .pdbqt and ligands must be .pdbqt", file=sys.stderr)
        print("  Install Meeko: pip install meeko", file=sys.stderr)
        print("  Or install Open Babel: apt install openbabel", file=sys.stderr)

    targets_to_run = list(requested_targets)
    if args.continue_run:
        pending_targets = [
            name for name in targets_to_run
            if not target_already_done(output_dir, name, method_label)
        ]
        skipped_count = len(targets_to_run) - len(pending_targets)
        if skipped_count:
            print(f"[continue] skipping {skipped_count} completed target(s) for {method_label}")
        targets_to_run = pending_targets

    if args.max_targets_per_run is not None:
        targets_to_run = targets_to_run[: args.max_targets_per_run]

    if not targets_to_run:
        refreshed = write_method_outputs(
            dataset_root=dataset_root,
            output_dir=output_dir,
            method_label=method_label,
            target_names=tracked_targets,
            mode=args.mode,
            center_source=args.center_source,
            seed=args.seed,
            invocation_wall_time_s=0.0,
        )
        print("No targets selected for execution after applying --continue/--max-targets-per-run.")
        print(f"Status:  {refreshed['paths']['status']}")
        print(f"Summary: {refreshed['paths']['summary']}")
        return refreshed["summary"]

    # ── run per target ──
    start_total = time.time()

    # subsampling workspace
    subsample_dir = output_dir / "_subsampled" if args.max_ligands else None

    for i, target_name in enumerate(targets_to_run, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(targets_to_run)}] TARGET: {target_name.upper()}")
        print(f"{'='*60}")
        target_result_dir = output_dir / target_name / method_label
        metrics_path = target_result_dir / "metrics.json"
        failure_path = target_result_dir / "failure.json"

        start_target = time.time()

        # determine effective dataset root for this target
        effective_root = dataset_root
        if args.max_ligands:
            print(f"  Subsampling to {args.max_ligands} ligands per label...")
            assert subsample_dir is not None
            prepare_subsampled_target(
                dataset_root, target_name, args.max_ligands, subsample_dir, seed=args.seed
            )
            effective_root = subsample_dir

        # build the dude_docking_benchmark.py command
        # NOTE: dude_docking_benchmark.py creates output_dir/target_name/ internally,
        # so we pass output_dir directly. Files will be moved post-run.
        bench_cmd = [
            sys.executable,
            str(BENCHMARKS_DIR / "dude_docking_benchmark.py"),
            "--dataset-root", str(effective_root),
            "--target", target_name,
            "--output-dir", str(output_dir),
            "--mode", args.mode,
            "--center-source", args.center_source,
            "--seed", str(args.seed),
            "--vina-cpu", str(args.vina_cpu),
            "--vina-exhaustiveness", str(args.vina_exhaustiveness),
            "--vina-num-modes", str(args.vina_num_modes),
            "--autosites", str(args.autosites),
            "--site-hit-threshold", str(args.site_hit_threshold),
        ]
        if rec_cmd:
            bench_cmd.extend(["--receptor-prepare-command", rec_cmd])
        if lig_cmd:
            bench_cmd.extend(["--ligand-prepare-command", lig_cmd])
        if args.bedroc_alpha is not None:
            bench_cmd.extend(["--bedroc-alpha", str(args.bedroc_alpha)])

        try:
            log_path = target_result_dir / "benchmark_runner.log"
            target_result_dir.mkdir(parents=True, exist_ok=True)
            with log_path.open("w", encoding="utf-8") as log_fh:
                proc = subprocess.run(
                    bench_cmd,
                    cwd=REPO_ROOT,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=args.timeout,
                )

            elapsed = time.time() - start_target

            if proc.returncode != 0:
                # read the log tail for a useful error message
                log_tail = ""
                try:
                    log_tail = log_path.read_text("utf-8")[-500:]
                except Exception:
                    pass
                print(f"  FAILED (exit code {proc.returncode}) in {elapsed:.0f}s")
                print(f"    See log: {log_path}")
                if log_tail:
                    print(f"    Last output: ...{log_tail.strip()[-200:]}")
                write_failure_record(
                    failure_path,
                    {
                        "target": target_name,
                        "method_label": method_label,
                        "status": "exit_code",
                        "reason": f"exit code {proc.returncode}",
                        "elapsed_s": round(elapsed, 3),
                        "log_path": str(log_path.resolve()),
                    },
                )
                refreshed = write_method_outputs(
                    dataset_root=dataset_root,
                    output_dir=output_dir,
                    method_label=method_label,
                    target_names=tracked_targets,
                    mode=args.mode,
                    center_source=args.center_source,
                    seed=args.seed,
                    invocation_wall_time_s=time.time() - start_total,
                )
                print(f"  Status updated: {refreshed['paths']['status']}")
                continue

            # read the summary.json produced by the subprocess
            summary_path = output_dir / "summary.json"
            target_entry = None
            sub_failed = False
            if summary_path.exists():
                summary = json.loads(summary_path.read_text("utf-8"))
                # check if this target is in the completed list
                for entry in summary.get("targets", []):
                    if entry.get("target") == target_name:
                        target_entry = entry
                        break
                # also check the failed list from dude_docking_benchmark.py
                for fail in summary.get("failed", []):
                    if fail.get("target") == target_name:
                        sub_failed = True
                        reason = fail.get("reason", "unknown error")
                        print(f"  FAILED in {elapsed:.0f}s")
                        # show a concise reason
                        reason_short = reason[:200] + "..." if len(reason) > 200 else reason
                        print(f"    Reason: {reason_short}")
                        print(f"    See log: {log_path}")
                        write_failure_record(
                            failure_path,
                            {
                                "target": target_name,
                                "method_label": method_label,
                                "status": "benchmark_failed",
                                "reason": reason_short,
                                "elapsed_s": round(elapsed, 3),
                                "log_path": str(log_path.resolve()),
                            },
                        )
                        break

            if sub_failed:
                refreshed = write_method_outputs(
                    dataset_root=dataset_root,
                    output_dir=output_dir,
                    method_label=method_label,
                    target_names=tracked_targets,
                    mode=args.mode,
                    center_source=args.center_source,
                    seed=args.seed,
                    invocation_wall_time_s=time.time() - start_total,
                )
                print(f"  Status updated: {refreshed['paths']['status']}")
                continue

            if target_entry:
                target_entry["runner_elapsed_s"] = round(elapsed, 3)
                target_entry["method_label"] = method_label
                target_entry["center_source"] = args.center_source
                target_entry["benchmark_runner_log"] = str(log_path.resolve())
                auc = target_entry.get("roc_auc", "?")
                auc_str = f"{auc:.4f}" if isinstance(auc, (int, float)) else "?"
                print(f"  [OK] DONE in {elapsed:.0f}s -- ROC-AUC = {auc_str}")
            else:
                # subprocess exited 0 but no results found — check the log
                log_tail = ""
                try:
                    log_tail = log_path.read_text("utf-8")[-500:]
                except Exception:
                    pass
                print(f"  FAILED in {elapsed:.0f}s — no results found")
                print(f"    See log: {log_path}")
                if log_tail:
                    print(f"    Last output: ...{log_tail.strip()[-200:]}")
                write_failure_record(
                    failure_path,
                    {
                        "target": target_name,
                        "method_label": method_label,
                        "status": "missing_summary_entry",
                        "reason": "no results in summary.json",
                        "elapsed_s": round(elapsed, 3),
                        "log_path": str(log_path.resolve()),
                    },
                )

            # ── Post-run: copy DB and move result files into method subfolder ──
            if target_entry:
                # Copy the SQLite results database into the method subfolder
                docking_db = REPO_ROOT / "docking" / "RESULTS_DIR" / "ultidock_results.db"
                if docking_db.exists():
                    dest_db = target_result_dir / "ultidock_results.db"
                    try:
                        copy_replace(docking_db, dest_db)
                        print(f"  DB copied to {dest_db}")
                    except Exception as e:
                        print(f"  DB copy failed: {e}")

                # Move result files from <target>/ into <target>/ultidock_<mode>/ if needed
                flat_target_dir = output_dir / target_name
                for fname in ("scores.csv", "setup.log", "dock_v02.log", "site_summary.json", "site_summary.csv", "site_metrics.json"):
                    src = flat_target_dir / fname
                    dst = target_result_dir / fname
                    if src.exists():
                        move_replace(src, dst)

                sites_src = flat_target_dir / "sites"
                sites_dst = target_result_dir / "sites"
                if sites_src.exists():
                    move_replace(sites_src, sites_dst)

                # Copy manifest files for later score recalculation
                docking_root = REPO_ROOT / "docking"
                for manifest in ("actives_manifest.json", "decoys_manifest.json"):
                    src = docking_root / manifest
                    dst = target_result_dir / manifest
                    if src.exists():
                        copy_replace(src, dst)

                # Copy auto-detected binding site centers for reproducibility
                centers_src = docking_root / "MACRO_MOL_DIR" / "receptor" / "centers.tsv"
                centers_dst = target_result_dir / "centers.tsv"
                if centers_src.exists():
                    copy_replace(centers_src, centers_dst)
                    print(f"  centers.tsv backed up to {centers_dst}")

                rewrite_metric_paths_for_method_dir(target_entry, target_result_dir)
                metrics_path.write_text(json.dumps(target_entry, indent=2), "utf-8")
                if failure_path.exists():
                    failure_path.unlink()

                # Generate per-target plots
                per_target_summary = target_result_dir / "summary.json"
                per_target_summary.write_text(json.dumps(
                    {"targets": [target_entry]}, indent=2
                ), "utf-8")

                try:
                    from benchmarks.dude_plots import generate_all_plots as _gap
                except ImportError:
                    try:
                        from dude_plots import generate_all_plots as _gap
                    except ImportError:
                        _gap = None
                if _gap:
                    try:
                        _gap(per_target_summary)
                        print(f"  Plots generated in {target_result_dir / 'plots'}")
                    except Exception as e:
                        print(f"  Plot generation failed: {e}")

            refreshed = write_method_outputs(
                dataset_root=dataset_root,
                output_dir=output_dir,
                method_label=method_label,
                target_names=tracked_targets,
                mode=args.mode,
                center_source=args.center_source,
                seed=args.seed,
                invocation_wall_time_s=time.time() - start_total,
            )
            print(f"  Status updated: {refreshed['paths']['status']}")

        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after {args.timeout}s")
            write_failure_record(
                failure_path,
                {
                    "target": target_name,
                    "method_label": method_label,
                    "status": "timeout",
                    "reason": "timeout",
                    "elapsed_s": round(time.time() - start_target, 3),
                    "log_path": str((target_result_dir / "benchmark_runner.log").resolve()),
                },
            )
            refreshed = write_method_outputs(
                dataset_root=dataset_root,
                output_dir=output_dir,
                method_label=method_label,
                target_names=tracked_targets,
                mode=args.mode,
                center_source=args.center_source,
                seed=args.seed,
                invocation_wall_time_s=time.time() - start_total,
            )
            print(f"  Status updated: {refreshed['paths']['status']}")
        except Exception as exc:
            print(f"  ERROR: {exc}")
            traceback.print_exc()
            write_failure_record(
                failure_path,
                {
                    "target": target_name,
                    "method_label": method_label,
                    "status": "runner_exception",
                    "reason": str(exc),
                    "elapsed_s": round(time.time() - start_target, 3),
                    "log_path": str((target_result_dir / "benchmark_runner.log").resolve()),
                },
            )
            refreshed = write_method_outputs(
                dataset_root=dataset_root,
                output_dir=output_dir,
                method_label=method_label,
                target_names=tracked_targets,
                mode=args.mode,
                center_source=args.center_source,
                seed=args.seed,
                invocation_wall_time_s=time.time() - start_total,
            )
            print(f"  Status updated: {refreshed['paths']['status']}")

    total_time = time.time() - start_total
    refreshed = write_method_outputs(
        dataset_root=dataset_root,
        output_dir=output_dir,
        method_label=method_label,
        target_names=tracked_targets,
        mode=args.mode,
        center_source=args.center_source,
        seed=args.seed,
        invocation_wall_time_s=total_time,
    )
    final_summary = refreshed["summary"]
    status = refreshed["status"]
    paths = refreshed["paths"]

    # ── aggregate results ──
    print(f"\n{'='*60}")
    print("BENCHMARK COMPLETE")
    print(f"{'='*60}")
    print(f"Total time: {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"Completed:  {status['n_completed']}/{status['n_targets']}")
    print(f"Failed:     {status['n_failed']}")
    print(f"Pending:    {status['n_pending']}")

    final_path = paths["summary"]
    compatibility_path = output_dir / "summary.json"
    compatibility_path.write_text(json.dumps(final_summary, indent=2), "utf-8")
    print(f"\nSummary: {final_path}")
    print(f"Status:  {paths['status']}")
    print(f"Runtime: {paths['runtime_summary']}")
    print(f"Errors:  {paths['failure_summary']}")

    # print mean metrics
    aggregate = final_summary.get("aggregate", {})
    if aggregate:
        print("\nAggregate metrics:")
        for k, v in aggregate.items():  # type: ignore[union-attr]
            print(f"  {k}: {v:.4f}")

    site_validation = final_summary.get("site_validation", {})
    if site_validation:
        print("\nSite validation:")
        for key in (
            "any_site_success_rate",
            "mean_best_distance_a",
            "mean_best_site_order",
        ):
            value = site_validation.get(key)
            if isinstance(value, float):
                print(f"  {key}: {value:.4f}")

    # ── generate plots ──
    try:
        from benchmarks.dude_plots import generate_all_plots
    except ImportError:
        try:
            from dude_plots import generate_all_plots
        except ImportError:
            generate_all_plots = None  # type: ignore

    if generate_all_plots and final_summary.get("targets"):
        print("\nGenerating scientific plots...")
        plots_dir = output_dir / f"{method_label}__plots"
        plots = generate_all_plots(final_path, plots_dir)
        if plots:
            print(f"  {len(plots)} plot(s) generated in {plots_dir}")

    return final_summary


# ── CLI ───────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the full DUD-E benchmark autonomously — download, dock, score, plot.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--targets",
        help="Comma-separated target names, or 'all'. Default: auto-discover from dataset root.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download missing DUD-E datasets before running.",
    )
    parser.add_argument(
        "--dataset-root",
        default=str(DEFAULT_DATASET_ROOT),
        help="Root directory containing target subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Directory where results, logs, and plots are written.",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "gpu", "cpu", "cuda", "opencl"],
        default="cpu",
        help="Docking mode.",
    )
    parser.add_argument(
        "--center-source",
        choices=["auto", "crystal"],
        default="auto",
        help="How the benchmark arm should obtain centers.tsv.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--vina-cpu", type=int, default=2, help="Vina CPU count.")
    parser.add_argument("--vina-exhaustiveness", type=int, default=8, help="Vina exhaustiveness.")
    parser.add_argument("--vina-num-modes", type=int, default=1, help="Vina num_modes.")
    parser.add_argument("--autosites", type=int, default=6, help="Auto-site count for center-source=auto.")
    parser.add_argument("--site-hit-threshold", type=float, default=4.0, help="Hit threshold in Angstrom for site evaluation.")
    parser.add_argument("--bedroc-alpha", type=float, default=20.0, help="BEDROC alpha parameter.")
    parser.add_argument(
        "--max-ligands",
        type=int,
        help="Subsample to N actives + N decoys per target (for quick tests).",
    )
    parser.add_argument(
        "--continue",
        dest="continue_run",
        action="store_true",
        help="Skip targets that already have scores.csv.",
    )
    parser.add_argument(
        "--max-targets-per-run",
        type=int,
        help="Run at most N selected targets this session. Useful for chunked receptor-by-receptor runs.",
    )
    parser.add_argument(
        "--refresh-summary",
        action="store_true",
        help="Rebuild the arm-level summary/status files from existing per-target outputs without docking.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate targets and configuration without docking.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Per-target timeout in seconds (default: no timeout).",
    )
    parser.add_argument(
        "--receptor-prepare-command",
        help="Override receptor conversion command template (e.g. for Meeko).",
    )
    parser.add_argument(
        "--ligand-prepare-command",
        help="Override ligand conversion command template (e.g. for Meeko).",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_benchmark_for_targets(args)


if __name__ == "__main__":
    main()
