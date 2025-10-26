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

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser used by both setup.py and run.py."""

    parser = argparse.ArgumentParser(description="Ultidock setup")
    parser.add_argument(
        "--wget",
        metavar="PATH",
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
        help="Choose the GPU detection mode (default: prompt)",
    )
    parser.add_argument(
        "--ligands-dir",
        "--LIGANDS_DIR",
        dest="ligands_dir",
        help="Directory to store ligand files",
    )
    parser.add_argument(
        "--docking-dir",
        "--DOCKING_DIR",
        dest="docking_dir",
        help="Directory to store docking outputs",
    )
    parser.add_argument(
        "--analysis-dir",
        "--ANALYSIS_DIR",
        dest="analysis_dir",
        help="Directory to store analysis artifacts",
    )
    parser.add_argument(
        "--vina-dir",
        "--VINA_DIR",
        dest="vina_dir",
        help="Directory containing AutoDock Vina binaries",
    )
    parser.add_argument(
        "--autodock-gpu-dir",
        "--AUTODOCK_GPU_DIR",
        dest="autodock_gpu_dir",
        help="Directory containing AutoDock-GPU",
    )
    parser.add_argument(
        "--macro-mol-dir",
        "--MACRO_MOL_DIR",
        dest="macro_mol_dir",
        help="Directory containing receptor PDBQT files",
    )
    parser.add_argument(
        "--results-dir",
        "--RESULTS_DIR",
        dest="results_dir",
        help="Directory to store docking results database",
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

    if not sys.stdin.isatty():
        print("[setup] Non-interactive session detected; skipping ligand download prompt.")
        return None

    prompt = "Enter the path to the .wget file"
    default = current_dir / "ligands.wget"
    return str(_normalize_path(ask_for_input(prompt, str(default))))

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
            print(f"[BUILD] Compiling AutoDock-GPU for {GPU_TYPE} (DEVICE={device_env}, NUMWI={NUMWI})…")
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
        print("[BUILD] CPU mode: building/checking AutoGrid only…")
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
    CURRENT_DIR = Path(__file__).resolve().parent
    print("=" * 50)
    print("Welcome to the Ultidock Setup")
    print("=" * 50)

    if args.mode:
        mode = _normalize_mode(args.mode)
    else:
        user_mode = input(
            "Select run mode (GPU = NVidia-CUDA/OpenCL, AMD-OpenCL or CPU or auto) [default: auto]: "
        )
        mode = _normalize_mode(user_mode or "auto")

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
        config_file.write('GRID_MODE = "centers"      # ligand | residues | centers | blind\n')
        config_file.write('GRID_SPACING = 0.375\n')
        config_file.write('GRID_MARGIN = 5.0         # Å\n')
        config_file.write('GRID_CAP = 30.0           # Å cap per axis for blind mode\n')
        config_file.write('AUTO_GRID_BIN = os.path.join(AUTODOCK_GPU_DIR, "autogrid", "autogrid4")\n')
        config_file.write('CENTERS_TSV  = os.path.join(MACRO_MOL_DIR, "centers.tsv")  # path or None\n')
        config_file.write('REF_LIGAND_PDB = None    # path to co-crystal/ref ligand if GRID_MODE="ligand"\n')
        config_file.write('HOTSPOT_NMS_MINSEP_A = 2.0\n')
        config_file.write('R_MIN_CAVITY_A = 20.0   # minimum inscribed sphere radius for cavity acceptance\n')
        config_file.write('SURFACE_SHELL__MIN_A = 2.0  # min/max distance from protein surface for surface pockets\n')
        config_file.write('SURFACE_SHELL__MAX_A = 20.0\n')
        config_file.write('SURFACE_NMS_MINSEP_A = 5        # voxels for non-max suppression of surface pockets\n')
        config_file.write('MAX_CENTER_DIST_A = 10.0         \n')
        config_file.write('CONTACT_SHELL_A = 4.0           # voxels ≤5 Å from surface count as “contact”\n')
        config_file.write('HOTSPOT_BOX_ANGLE = 35       # minimum box side length (Å)\n')
        config_file.write('MIN_SURFACE_FRAC = 0.01        # ~0.2% of box must be near-surface\n')
        config_file.write('AUTOSITES = 6\n')
        print(f"Configuration saved to {config_path} and directories were ensured!")

    if wget_file_path:
        download_ligands_from_file(wget_file_path, ligands_dir)

    return config_values


def main(argv: Optional[Iterable[str]] = None) -> dict:
    args = parse_args(argv)
    return run_setup(args)


if __name__ == "__main__":
    main()
