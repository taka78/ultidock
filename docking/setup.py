import os
import subprocess
import sys
import shutil
import stat
import argparse
from pathlib import Path
from typing import Iterable, Optional

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
NUMWI = "128"  # Default number of work items
CURRENT_DIR = Path(SCRIPT_DIR)
DEFAULT_LIGANDS_DIR = str((CURRENT_DIR / "LIGANDS_DIR").resolve())
DEFAULT_DOCKING_DIR = str((CURRENT_DIR / "DOCKING_DIR").resolve())
DEFAULT_ANALYSIS_DIR = str((CURRENT_DIR / "ANALYSIS_DIR").resolve())
DEFAULT_VINA_DIR = str((CURRENT_DIR / "VINA_DIR").resolve())
DEFAULT_AUTODOCK_GPU_DIR = str((CURRENT_DIR / "AUTODOCK_GPU_DIR").resolve())
DEFAULT_MACRO_MOL_DIR = str((CURRENT_DIR / "MACRO_MOL_DIR").resolve())
DEFAULT_RESULTS_DIR = str((CURRENT_DIR / "RESULTS_DIR").resolve())
DEFAULT_CENTERS_TSV = str((CURRENT_DIR / "MACRO_MOL_DIR" / "centers.tsv").resolve())
DEFAULT_WGET_PATH = str((CURRENT_DIR / "ligands.wget").resolve())

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser used by both setup.py and run.py."""

    parser = argparse.ArgumentParser(description="Ultidock setup")
    parser.add_argument(
        "--wget",
        metavar="PATH",
        default=DEFAULT_WGET_PATH,
        help="Path to a .wget file (overrides prompt)",
    )
    parser.add_argument(
        "--skip-wget",
        action="store_true",
        help="Skip executing wget commands (ligands must already be present)",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "gpu", "cpu", "cuda", "opencl"],
        default="auto",
        help="Choose the GPU detection mode",
    )
    parser.add_argument(
        "--ligands-dir",
        "--LIGANDS_DIR",
        dest="ligands_dir",
        default=DEFAULT_LIGANDS_DIR,
        help="Directory to store ligand files",
    )
    parser.add_argument(
        "--docking-dir",
        "--DOCKING_DIR",
        dest="docking_dir",
        default=DEFAULT_DOCKING_DIR,
        help="Directory to store docking outputs",
    )
    parser.add_argument(
        "--analysis-dir",
        "--ANALYSIS_DIR",
        dest="analysis_dir",
        default=DEFAULT_ANALYSIS_DIR,
        help="Directory to store analysis artifacts",
    )
    parser.add_argument(
        "--vina-dir",
        "--VINA_DIR",
        dest="vina_dir",
        default=DEFAULT_VINA_DIR,
        help="Directory containing AutoDock Vina binaries",
    )
    parser.add_argument(
        "--autodock-gpu-dir",
        "--AUTODOCK_GPU_DIR",
        dest="autodock_gpu_dir",
        default=DEFAULT_AUTODOCK_GPU_DIR,
        help="Directory containing AutoDock-GPU",
    )
    parser.add_argument(
        "--macro-mol-dir",
        "--MACRO_MOL_DIR",
        dest="macro_mol_dir",
        default=DEFAULT_MACRO_MOL_DIR,
        help="Directory containing receptor PDBQT files",
    )
    parser.add_argument(
        "--results-dir",
        "--RESULTS_DIR",
        dest="results_dir",
        default=DEFAULT_RESULTS_DIR,
        help="Directory to store docking results database",
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Benchmark mode: skip wget prompt unless explicitly provided and persist benchmark-oriented config values.",
    )
    parser.add_argument(
        "--grid-mode",
        choices=["ligand", "residues", "centers", "blind"],
        default="centers",
        help="Grid generation mode written to config.py.",
    )
    parser.add_argument(
        "--grid-spacing",
        type=float,
        default=0.375,
        help="Grid spacing written to config.py.",
    )
    parser.add_argument(
        "--grid-margin",
        type=float,
        default=5.0,
        help="Grid margin written to config.py.",
    )
    parser.add_argument(
        "--grid-cap",
        type=float,
        default=150.0,
        help="Grid cap written to config.py.",
    )
    parser.add_argument(
        "--centers-tsv",
        default=DEFAULT_CENTERS_TSV,
        help="Path to a centers.tsv file written to config.py.",
    )
    parser.add_argument(
        "--ref-ligand-pdb",
        help='Reference ligand path for GRID_MODE="ligand".',
    )
    parser.add_argument(
        "--autosites",
        type=int,
        default=6,
        help="Number of automatically generated sites written to config.py.",
    )
    parser.add_argument(
        "--vina-cpu",
        type=int,
        default=2,
        help="CPU count passed to Vina in CPU mode.",
    )
    parser.add_argument(
        "--vina-seed",
        type=int,
        help="Seed passed to Vina in CPU mode.",
    )
    parser.add_argument(
        "--vina-exhaustiveness",
        type=int,
        default=8,
        help="Exhaustiveness passed to Vina in CPU mode.",
    )
    parser.add_argument(
        "--vina-num-modes",
        type=int,
        default=9,
        help="num_modes passed to Vina in CPU mode.",
    )
    parser.add_argument(
        "--skip-profile",
        action="store_true",
        help="Skip per-receptor profiling (preserves hand-tuned .config.toml sidecars).",
    )
    parser.add_argument(
        "--receptor-prep-mode",
        choices=["auto", "off"],
        default="auto",
        help="Prepare/canonicalize receptor inputs in MACRO_MOL_DIR using molguard.",
    )
    parser.add_argument(
        "--receptor-prepare-command",
        help="Override receptor conversion command template; use {input}, {output}, {seed}.",
    )
    parser.add_argument(
        "--receptor-prep-seed",
        type=int,
        default=42,
        help="Seed forwarded to receptor-prep conversion commands.",
    )
    parser.add_argument(
        "--force-receptor-prep",
        action="store_true",
        help="Regenerate receptor PDBQT outputs even when sibling outputs already exist.",
    )
    return parser


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments, optionally from an iterable."""

    parser = build_parser()
    return parser.parse_args(argv)



def ask_for_input(prompt, default):
    user_input = input(f"{prompt} [default: {default}]: ")
    return user_input if user_input else default


def _normalize_path(value: str) -> Path:
    """Expand user/home markers and resolve the provided path."""

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path.resolve()


def create_directory_if_needed(path: str | Path) -> str:
    """Ensure a directory exists exactly where the user requested it."""

    directory = Path(path).expanduser()
    if not directory.is_absolute():
        directory = Path.cwd() / directory
    directory.mkdir(parents=True, exist_ok=True)
    print(f"[ok] ensured directory {directory}")
    return str(directory)


def check_and_fix_receptors(
    macro_mol_dir: str | Path,
    *,
    prep_mode: str = "auto",
    prepare_command: str | None = None,
    seed: int = 42,
    force: bool = False,
) -> None:
    """
    Prepare receptor inputs through molguard, then scan every .pdbqt file,
    print a grouped summary of issues, then ask the user whether to
    auto-fix the problematic files with canonicalize_receptor.

    Requires molguard (pip install -e . from repo root).
    Silently skips if molguard is not importable (e.g. bare environment).
    """
    try:
        from molguard.io.pdbqt import LintError, canonicalize_receptor, pdbqt_check
        from molguard.io.receptor_prep import prepare_receptors_in_directory
    except ImportError:
        print("[warn] molguard not found -- skipping receptor prep/PDBQT check.")
        print("       Install with: pip install -e . (from repo root)")
        return

    mol_dir = Path(macro_mol_dir)
    if prep_mode.lower() not in {"off", "skip", "none", "false", "0"}:
        records = prepare_receptors_in_directory(
            mol_dir,
            prepare_command=prepare_command,
            seed=seed,
            timestamp="SETUP",
            force=force,
        )
        if records:
            print(f"[setup] molguard prepared/canonicalized {len(records)} receptor input(s).")
            for record in records:
                print(
                    f"  [receptor-prep] {record.output_path.name}  {record.action}  "
                    f"sha256={record.digest[:12]}..."
                )
    else:
        print("[setup] receptor prep disabled by --receptor-prep-mode=off.")

    pdbqt_files = sorted(mol_dir.rglob("*.pdbqt"))

    if not pdbqt_files:
        print(f"[info] No PDBQT files found in {mol_dir} -- nothing to check.")
        return

    print(f"\n=== Receptor PDBQT check ({len(pdbqt_files)} file(s) in {mol_dir}) ===")

    # Collect results: list of (path, report)
    issues: list[tuple[Path, object]] = []
    clean: list[Path] = []

    for pdbqt in pdbqt_files:
        report = pdbqt_check(pdbqt)
        if report.ok and not report.warnings:
            clean.append(pdbqt)
        else:
            issues.append((pdbqt, report))

    # Summary header
    print(f"  Clean  : {len(clean)} file(s)")
    print(f"  Issues : {len(issues)} file(s)")

    if not issues:
        print("  All receptor files passed the PDBQT check.")
        return

    # Print per-file issue details
    print()
    for pdbqt, report in issues:
        rel = pdbqt.relative_to(mol_dir) if pdbqt.is_relative_to(mol_dir) else pdbqt.name
        print(f"  [{rel}]")
        for e in report.errors:
            print(f"    ERROR   line {e.line_no:4d} [{e.column}] {e.code}: {e.message}")
        for w in report.warnings:
            print(f"    WARNING line {w.line_no:4d} [{w.column}] {w.code}: {w.message}")
        if report.errors:
            print(f"    -> {len(report.errors)} error(s), {len(report.warnings)} warning(s)")
        else:
            print(f"    -> 0 errors, {len(report.warnings)} warning(s) only")
    print()

    # Ask the user whether to auto-fix the problematic files
    try:
        answer = input(
            f"Auto-fix {len(issues)} file(s) with receptor canonicalization? "
            "This rewrites numeric columns and sorts atoms deterministically.\n"
            "WARNING: the original files will be overwritten in-place.\n"
            "[y/N]: "
        ).strip().lower()
    except (EOFError, KeyboardInterrupt):
        # Non-interactive environment (CI, piped stdin) -- skip without crashing
        print("[info] Non-interactive mode -- skipping auto-fix.")
        return

    if answer not in ("y", "yes"):
        print("[info] Skipped. You can fix individual files later with:")
        print("       ultidock pdbqt canonicalize-receptor <file> -o <file>")
        return

    fixed = 0
    failed = 0
    for pdbqt, report in issues:
        # Canonicalize if there are errors OR any BFAC_OVERFLOW warnings.
        # BFAC_OVERFLOW is a warning (not an error) because the bytes are valid,
        # but autogrid4's whitespace-token parser still misparsed the line —
        # so we must reformat even when there are no other structural errors.
        has_bfac_overflow = any(w.code == "BFAC_OVERFLOW" for w in report.warnings)
        if not report.errors and not has_bfac_overflow:
            print(f"  [skip] {pdbqt.name} -- warnings only, no fix needed")
            continue
        reason = "errors" if report.errors else "BFAC_OVERFLOW"
        try:
            digest = canonicalize_receptor(pdbqt, pdbqt, timestamp="SETUP")
            print(f"  [fixed] {pdbqt.name}  ({reason})  sha256={digest[:12]}...")
            fixed += 1
        except (LintError, Exception) as exc:
            print(f"  [FAIL]  {pdbqt.name}  could not fix: {exc}")
            failed += 1

    print(f"\n  Fixed: {fixed}  Failed: {failed}")
    if failed:
        print("  [warn] Some files could not be fixed automatically.")
        print("         Check them with: ultidock pdbqt check <file>")


def detect_gpu():
    """Detect the GPU type (NVIDIA, AMD, or CPU fallback)."""
    
    # Check for NVIDIA GPU using nvidia-smi
    has_nvidia = shutil.which("nvidia-smi")
    has_amd = shutil.which("rocm-smi")

    if has_nvidia:
        try:
            subprocess.run(["nvidia-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            print("NVIDIA GPU detected.")
            return "CUDA"
        except subprocess.CalledProcessError:
            print("NVIDIA GPU detected but inaccessible.")

    if has_amd:
        try:
            subprocess.run(["rocm-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            print("AMD GPU detected.")
            return "OPENCL"
        except subprocess.CalledProcessError:
            print("AMD GPU detected but inaccessible.")

    print("No compatible GPU detected. Using CPU mode.")
    return "CPU"

def _normalize_mode(s: str) -> str:
    """
    Return one of: 'auto' | 'gpu' | 'cpu' | 'opencl' | 'cuda'
    Matches keywords anywhere in the input (case-insensitive).
    Priority: last matching keyword wins if multiple appear.
    """
    t = (s or "").lower()

    # (keyword -> normalized mode)
    keywords = [
        # explicit backends
        ("opencl", "opencl"),
        ("opengl", "opencl"),   # common slip
        ("ocl",    "opencl"),
        ("rocm",   "opencl"),
        ("amd",    "opencl"),

        ("cuda",   "cuda"),
        ("nvidia", "cuda"),
        (" nv ",   "cuda"),     # crude guard to avoid matching 'env'
        (" cu ",   "cuda"),

        # generic modes
        ("cpu",    "cpu"),
        ("gpu",    "gpu"),
        ("auto",   "auto"),
    ]

    last_hit = None
    last_pos = -1
    tt = f" {t} "  # pad to make the ' nv ' / ' cu ' checks work
    for kw, mode in keywords:
        pos = tt.rfind(kw)  # last occurrence wins
        if pos != -1 and pos > last_pos:
            last_pos = pos
            last_hit = mode

    return last_hit or "auto"



def download_ligands_from_file(wget_file_path, LIGANDS_DIR):
    if not os.path.exists(wget_file_path):
        print(f"Error: {wget_file_path} does not exist.")
        return

    print(f"Executing wget commands from {wget_file_path} to {LIGANDS_DIR}...")

    try:
        with open(wget_file_path, 'r') as file:
            for command in file:
                command = command.strip()
                if command:
                    if '-O' in command:
                        parts = command.split('-O')
                        url = parts[0].strip()
                        filename = parts[1].strip()
                        full_output_path = os.path.join(LIGANDS_DIR, filename)
                        final_command = f"{url} -O {full_output_path}"
                    else:
                        final_command = f"{command} -P {LIGANDS_DIR}"

                    print(f"Executing: {final_command}")
                    subprocess.run(final_command, shell=True, check=True)

        print(f"All downloads completed and stored in {LIGANDS_DIR}")
    except subprocess.CalledProcessError as e:
        print(f"Error during command execution: {e}")


def _bin_backend(path: str) -> str | None:
    """Return 'CUDA' if the binary links to libcuda/libcudart, 'OPENCL' if it links to libOpenCL, else None."""
    try:
        out = subprocess.run(["ldd", path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=2)
        s = out.stdout or ""
    except Exception:
        return None
    if "libOpenCL.so" in s:
        return "OPENCL"
    if "libcuda.so" in s or "libcudart.so" in s:
        return "CUDA"
    return None


def _resolve_directory(prompt: str, default_path: Path, cli_value: Optional[str]) -> str:
    """Resolve a directory path, preferring CLI values over interactive input."""

    if cli_value:
        target = _normalize_path(cli_value)
    else:
        user_value = ask_for_input(prompt, str(default_path))
        target = _normalize_path(user_value)
    return create_directory_if_needed(target)


def _resolve_wget_path(current_dir: Path, args: argparse.Namespace) -> Optional[str]:
    if getattr(args, "skip_wget", False):
        print("[setup] --skip-wget specified; skipping ligand download step.")
        return None

    if args.wget:
        return str(_normalize_path(args.wget))

    if getattr(args, "benchmark", False):
        print("[setup] Benchmark mode: skipping ligand download prompt.")
        return None

    if not sys.stdin.isatty():
        print("[setup] Non-interactive session detected; skipping ligand download prompt.")
        return None

    return str(_normalize_path(DEFAULT_WGET_PATH))

def _find_autodock_gpu_bin(autodock_dir: str, gpu_type: str, numwi) -> str | None:
    """
    Find a usable AutoDock-GPU binary for the requested backend ('CUDA'|'OPENCL').
    Accept both generic name (autodock_gpu_<N>wi) and backend-specific names.
    Verify by checking linked libraries with ldd.
    """
    bin_dir = os.path.join(autodock_dir, "bin")
    numwi = str(numwi).lower()
    generic = f"autodock_gpu_{numwi}wi"

    # Try the common names first (including the generic one you’re using for OCL too)
    names = [
        generic,
        "autodock_gpu",                # some builds are multi-backend
        "autodock_gpu_cuda",
        "autodock_gpu_ocl",
        "autodock_gpu_opencl",
        f"autodock_gpu_{numwi}wi_ocl",
        f"autodock_gpu_ocl_{numwi}wi",
    ]
    for name in names:
        p = os.path.join(bin_dir, name)
        if not (os.path.isfile(p) and os.access(p, os.X_OK)):
            continue
        be = _bin_backend(p)
        if be is None:
            # Unknown, but if the name matches the generic target, accept it
            if os.path.basename(p) == generic:
                return p
            continue
        if be == gpu_type.upper():
            return p
    return None

def detect_and_compile_autodock_gpu(AUTODOCK_GPU_DIR, GPU_TYPE, NUMWI):
    """
    Always ensure AutoGrid exists. If GPU mode (CUDA/OPENCL), also ensure AutoDock-GPU.
    In CPU mode we skip GPU binary checks/build and just build/check AutoGrid.
    """
    gpu_upper = (GPU_TYPE or "CPU").upper()

    compiler_script = os.path.join(SCRIPT_DIR, "autodock-gpu-compiler.sh")
    if not os.path.isfile(compiler_script):
        print(f"Compiler script not found at {compiler_script}")
        sys.exit(1)

    # Map to script's DEVICE
    if gpu_upper in ("NVIDIA", "CUDA"):
        device_env = "CUDA"
    elif gpu_upper in ("OPENCL", "AMD"):
        device_env = "OPENCL"
    else:
        device_env = "CPU"

    env = os.environ.copy()
    env["DEVICE"] = device_env
    env["NUMWI"] = str(NUMWI)
    env["GPU_TYPE"] = GPU_TYPE or device_env
    env["GPU_BACKEND"] = device_env

    found_bin = None  # <- important: define upfront so we never reference an unbound name

    # If GPU mode, try to reuse or build AutoDock-GPU
    if device_env in ("CUDA", "OPENCL"):
        found_bin = _find_autodock_gpu_bin(AUTODOCK_GPU_DIR, device_env, NUMWI)
        if found_bin:
            print(f"AutoDock-GPU is set up for {GPU_TYPE}: {found_bin}")
        else:
            print(f"[BUILD] Compiling AutoDock-GPU for {GPU_TYPE} (DEVICE={device_env}, NUMWI={NUMWI})...")
            subprocess.run(
                [
                    "bash",
                    compiler_script,
                    AUTODOCK_GPU_DIR,
                    device_env,
                    str(NUMWI),
                    GPU_TYPE or device_env,
                ],
                check=True,
                cwd=AUTODOCK_GPU_DIR,
                env=env,
            )
            # Re-check after build using the finder so we accept legacy names too
            found_bin = _find_autodock_gpu_bin(AUTODOCK_GPU_DIR, device_env, NUMWI)
            if not found_bin:
                print("AutoDock-GPU compilation finished but no binary was found in bin/.")
                sys.exit(1)
            st = os.stat(found_bin)
            os.chmod(found_bin, st.st_mode | stat.S_IXUSR)
            print(f"[BUILD] AutoDock-GPU ready: {found_bin}")
    else:
        # CPU-only: still run the script so it compiles/checks AutoGrid
        print("[BUILD] CPU mode: building/checking AutoGrid only...")
        subprocess.run(
            ["bash", compiler_script, AUTODOCK_GPU_DIR],
            check=True, cwd=AUTODOCK_GPU_DIR, env=env
        )

    # In all modes, verify AutoGrid
    autogrid_bin = os.path.join(AUTODOCK_GPU_DIR, "autogrid", "autogrid4")
    if not (os.path.isfile(autogrid_bin) and os.access(autogrid_bin, os.X_OK)):
        raise RuntimeError(f"AutoGrid not found or not executable at {autogrid_bin}")
    print(f"[BUILD] AutoGrid ready: {autogrid_bin}")

    # Optional: return paths if you want to use them later
    return {"autogrid": autogrid_bin, "adgpu": found_bin}


'''    # Try compiling AutoDock-GPU if needed
    if GPU_TYPE in ("NVIDIA", "AMD"):
        AUTODOCK_GPU_BIN = os.path.join(AUTODOCK_GPU_DIR, "bin", "autodock_gpu_128wi")

        if not os.path.isfile(AUTODOCK_GPU_BIN) or not os.access(AUTODOCK_GPU_BIN, os.X_OK):
            print("AutoDock-GPU binary not found. Attempting to compile it...")
            compiler_script = os.path.join(SCRIPT_DIR, "autodock-gpu-compiler.sh")
            try:
                subprocess.run(["bash", compiler_script, AUTODOCK_GPU_DIR], check=True)
                print("AutoDock-GPU compilation successful.")
            except subprocess.CalledProcessError as e:
                print("AutoDock-GPU compilation failed.")
                print(e)
                sys.exit(1)
        else:
            print("AutoDock-GPU binary already compiled and executable.")'''


def run_setup(args: argparse.Namespace) -> dict:
    print("=" * 50)
    print("Welcome to the Ultidock Setup")
    print("=" * 50)

    mode = _normalize_mode(args.mode)

    if mode in ("opencl", "cuda", "cpu"):
        GPU_TYPE = mode.upper()
    else:  # 'gpu' or 'auto'
        GPU_TYPE = detect_gpu()

    ligands_dir = _resolve_directory(
        "Enter the path for ligand files",
        CURRENT_DIR / "LIGANDS_DIR",
        getattr(args, "ligands_dir", None),
    )
    docking_dir = _resolve_directory(
        "Enter the path for docking files",
        CURRENT_DIR / "DOCKING_DIR",
        getattr(args, "docking_dir", None),
    )
    analysis_dir = _resolve_directory(
        "Enter the path for analysis files",
        CURRENT_DIR / "ANALYSIS_DIR",
        getattr(args, "analysis_dir", None),
    )
    vina_dir = _resolve_directory(
        "Enter the path for Autodock Vina",
        CURRENT_DIR / "VINA_DIR",
        getattr(args, "vina_dir", None),
    )
    autodock_gpu_dir = _resolve_directory(
        "Enter the path for Autodock-GPU files",
        CURRENT_DIR / "AUTODOCK_GPU_DIR",
        getattr(args, "autodock_gpu_dir", None),
    )
    macro_mol_dir = _resolve_directory(
        "Enter the path for macro molecule of your choice",
        CURRENT_DIR / "MACRO_MOL_DIR",
        getattr(args, "macro_mol_dir", None),
    )
    results_dir = _resolve_directory(
        "Enter the path for results files",
        CURRENT_DIR / "RESULTS_DIR",
        getattr(args, "results_dir", None),
    )

    wget_file_path = _resolve_wget_path(CURRENT_DIR, args)

    print(f"Setting up AutoDock-GPU and/or Autogrid in {autodock_gpu_dir} for {GPU_TYPE} mode...")
    detect_and_compile_autodock_gpu(autodock_gpu_dir, GPU_TYPE, NUMWI)

    config_path = Path(ROOT_DIR) / "docking" / "config.py"
    config_values = {
        "LIGANDS_DIR": ligands_dir,
        "DOCKING_DIR": docking_dir,
        "ANALYSIS_DIR": analysis_dir,
        "VINA_DIR": vina_dir,
        "AUTODOCK_GPU_DIR": autodock_gpu_dir,
        "MACRO_MOL_DIR": macro_mol_dir,
        "RESULTS_DIR": results_dir,
        "GPU_TYPE": GPU_TYPE,
        "DB_PATH": str(Path(results_dir) / "ultidock_results.db"),
        "NUMWI": NUMWI,
        "GRID_MODE": args.grid_mode,
        "GRID_SPACING": args.grid_spacing,
        "GRID_MARGIN": args.grid_margin,
        "GRID_CAP": args.grid_cap,
        "CENTERS_TSV": (
            str(_normalize_path(args.centers_tsv))
            if args.centers_tsv
            else str(Path(macro_mol_dir) / "centers.tsv")
        ),
        "REF_LIGAND_PDB": (
            str(_normalize_path(args.ref_ligand_pdb))
            if args.ref_ligand_pdb
            else None
        ),
        "AUTOSITES": args.autosites,
        "VINA_CPU": args.vina_cpu,
        "VINA_SEED": args.vina_seed,
        "VINA_EXHAUSTIVENESS": args.vina_exhaustiveness,
        "VINA_NUM_MODES": args.vina_num_modes,
        "BENCHMARK_MODE": bool(args.benchmark),
        "RECEPTOR_PREP_MODE": args.receptor_prep_mode,
        "RECEPTOR_PREP_COMMAND": args.receptor_prepare_command,
        "RECEPTOR_PREP_SEED": args.receptor_prep_seed,
        "RECEPTOR_PREP_FORCE": bool(args.force_receptor_prep),
    }

    with open(config_path, "w", encoding="utf-8") as config_file:
        config_file.write("# config.py\n")
        config_file.write("# Auto-generated config.py\n")
        config_file.write("import os\n\n")
        config_file.write("BASE_DIR = os.path.abspath(os.path.dirname(__file__))\n\n")
        for key in (
            "LIGANDS_DIR",
            "DOCKING_DIR",
            "ANALYSIS_DIR",
            "VINA_DIR",
            "AUTODOCK_GPU_DIR",
            "MACRO_MOL_DIR",
            "RESULTS_DIR",
        ):
            config_file.write(f"{key} = {repr(config_values[key])}\n")
        config_file.write(f"GPU_TYPE = {repr(config_values['GPU_TYPE'])}\n")
        config_file.write(f"DB_PATH = {repr(config_values['DB_PATH'])}\n")
        config_file.write(f"NUMWI = {repr(config_values['NUMWI'])}\n")
        config_file.write(
            f"GRID_MODE = {repr(config_values['GRID_MODE'])}      # ligand | residues | centers | blind\n"
        )
        config_file.write("SITE_POLICY = 'receptor_search'  # receptor_search | exhaustive_search | internal | surface | hybrid\n")
        config_file.write(f"GRID_SPACING = {repr(config_values['GRID_SPACING'])}\n")
        config_file.write(f"GRID_MARGIN = {repr(config_values['GRID_MARGIN'])}         # Å\n")
        config_file.write(f"GRID_CAP = {repr(config_values['GRID_CAP'])}           # Å cap per axis for blind mode\n")
        config_file.write('AUTO_GRID_BIN = os.path.join(AUTODOCK_GPU_DIR, "autogrid", "autogrid4")\n')
        config_file.write(f"CENTERS_TSV = {repr(config_values['CENTERS_TSV'])}  # path or None\n")
        config_file.write(
            f"REF_LIGAND_PDB = {repr(config_values['REF_LIGAND_PDB'])}    # path to co-crystal/ref ligand if GRID_MODE='ligand'\n"
        )
        config_file.write('HOTSPOT_NMS_MINSEP_A = 14.0  # minimum center separation; profile_receptors may tune per receptor\n')
        config_file.write('R_MIN_CAVITY_A = None    # None = adaptive per receptor; set a float to force a fixed threshold\n')
        config_file.write('ADAPTIVE_R_MIN_PERCENTILE = 85.0\n')
        config_file.write('ADAPTIVE_R_MIN_PEAK_PERCENTILE = 10.0\n')
        config_file.write('ADAPTIVE_R_MIN_FLOOR_A = 2.0\n')
        config_file.write('ADAPTIVE_R_MIN_CEIL_A = 5.0\n')
        config_file.write('ADAPTIVE_R_MIN_POCKET_ZONE_MAX_A = 8.0\n')
        config_file.write('ADAPTIVE_R_MIN_PEAK_WINDOW_A = 4.0\n')
        config_file.write('MAPS_POCKET_MAX_A = 15.0\n')
        config_file.write('HOTSPOT_NMS_BOX_FRACTION = 0.40\n')
        config_file.write('HOTSPOT_NMS_MIN_A = 4.0\n')
        config_file.write('HOTSPOT_NMS_MAX_A = 25.0\n')
        config_file.write('SURFACE_SHELL__MIN_A = 2.0  # min/max distance from protein surface for surface pockets\n')
        config_file.write('SURFACE_SHELL__MAX_A = 20.0\n')
        config_file.write('SURFACE_NMS_MINSEP_A = 5.0        # Angstrom floor for surface-pocket non-max suppression\n')
        config_file.write('MAX_CENTER_DIST_A = 10.0         \n')
        config_file.write('CONTACT_SHELL_A = 4.0           # Angstrom shell near protein surface counted as contact\n')
        config_file.write('HOTSPOT_BOX_ANGLE = 35       # minimum box side length (Å)\n')
        config_file.write('MIN_SURFACE_FRAC = 0.01        # minimum local near-surface fraction for accepting surface pockets\n')
        config_file.write(f"AUTOSITES = {repr(config_values['AUTOSITES'])}\n")
        config_file.write(f"VINA_CPU = {repr(config_values['VINA_CPU'])}\n")
        config_file.write(f"VINA_SEED = {repr(config_values['VINA_SEED'])}\n")
        config_file.write(f"VINA_EXHAUSTIVENESS = {repr(config_values['VINA_EXHAUSTIVENESS'])}\n")
        config_file.write(f"VINA_NUM_MODES = {repr(config_values['VINA_NUM_MODES'])}\n")
        config_file.write(f"RECEPTOR_PREP_MODE = {repr(config_values['RECEPTOR_PREP_MODE'])}\n")
        config_file.write(f"RECEPTOR_PREP_COMMAND = {repr(config_values['RECEPTOR_PREP_COMMAND'])}\n")
        config_file.write(f"RECEPTOR_PREP_SEED = {repr(config_values['RECEPTOR_PREP_SEED'])}\n")
        config_file.write(f"RECEPTOR_PREP_FORCE = {repr(config_values['RECEPTOR_PREP_FORCE'])}\n")
        config_file.write(f"BENCHMARK_MODE = {repr(config_values['BENCHMARK_MODE'])}\n")
        print(f"Configuration saved to {config_path} and directories were ensured!")

    if wget_file_path:
        download_ligands_from_file(wget_file_path, ligands_dir)

    # Check receptor files in MACRO_MOL_DIR for AutoDock column-format issues.
    # This runs after the config is written so macro_mol_dir is confirmed.
    check_and_fix_receptors(
        macro_mol_dir,
        prep_mode=args.receptor_prep_mode,
        prepare_command=args.receptor_prepare_command,
        seed=args.receptor_prep_seed,
        force=args.force_receptor_prep,
    )

    # Auto-generate per-receptor .config.toml sidecars (unless suppressed).
    # --skip-profile preserves hand-tuned sidecars from previous runs.
    if not getattr(args, "skip_profile", False):
        try:
            from profile_receptors import profile_all
            profile_all(Path(macro_mol_dir), force=False, dry_run=False)
        except ImportError:
            print("[warn] profile_receptors.py not found — skipping per-receptor profiling.")
        except Exception as exc:
            print(f"[warn] Per-receptor profiling failed: {exc} — continuing without sidecars.")
    else:
        print("[setup] --skip-profile set — skipping per-receptor config generation.")

    return config_values


def main(argv: Optional[Iterable[str]] = None) -> dict:
    args = parse_args(argv)
    return run_setup(args)


if __name__ == "__main__":
    main()
