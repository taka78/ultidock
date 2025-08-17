import os
import subprocess
import sys
import shutil
import stat
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
NUMWI = "64"  # Default number of work items

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))


def ask_for_input(prompt, default):
    user_input = input(f"{prompt} [default: {default}]: ")
    return user_input if user_input else default


def create_directory_if_needed(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory {directory} created.")
    else:
        print(f"Directory {directory} already exists, continuing without creating it.")


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
    Ensures the correct AutoDock-GPU binary exists.
    - For NVIDIA: expects CUDA build (e.g., autodock_gpu_128wi).
    - For AMD/OPENCL: expects an OpenCL build (autodock_gpu_ocl* or autodock_gpu_cuda*).
    Calls ./autodock-gpu-compiler.sh to build if missing.
    """
    gpu_upper = (GPU_TYPE or "CPU").upper()
    needs_gpu = gpu_upper in ("CUDA", "OPENCL")
    if not needs_gpu:
        print("AutoDock will run in CPU mode.")
        return

    # 1) try to find an existing binary
    found = _find_autodock_gpu_bin(AUTODOCK_GPU_DIR, gpu_upper, NUMWI)
    if found:
        print(f"AutoDock-GPU is set up for {GPU_TYPE}: {found}")
        return

    # 2) build via your compiler script
    compiler_script = os.path.join(SCRIPT_DIR, "autodock-gpu-compiler.sh")
    if not os.path.isfile(compiler_script):
        print(f"Compiler script not found at {compiler_script}")
        sys.exit(1)

    # pick DEVICE for the script
    if gpu_upper == "NVIDIA" or gpu_upper == "CUDA":
        device_env = "CUDA"
    elif gpu_upper == "OPENCL" or gpu_upper == "AMD":
        device_env = "OPENCL"
    else:
        print(f"Unsupported GPU type: {GPU_TYPE}. Only CUDA and OPENCL are supported, cpu mode will be used.")
        device_env = "CPU"

    env = os.environ.copy()
    env["DEVICE"] = device_env
    env["NUMWI"] = str(NUMWI)

    print(f"[BUILD] AutoDock-GPU binary missing. Compiling for {GPU_TYPE} (DEVICE={device_env}, NUMWI={NUMWI})…")
    try:
        # Pass AUTODOCK_GPU_DIR as the script argument (your script expects it)
        print(env)
        subprocess.run(
            ["bash", compiler_script, AUTODOCK_GPU_DIR],
            check=True,
            cwd=AUTODOCK_GPU_DIR,
            env=env,
        )
    except subprocess.CalledProcessError:
        print("Compiler script failed.")
        sys.exit(1)

    # 3) re-check after build
    if device_env == "CUDA":
        found = Path(f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_cuda_{NUMWI}wi")
    elif device_env == "OPENCL":
        found = Path(f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_ocl_{NUMWI}wi")

    if found.exists():
        print(f"Found AutoDock-GPU binary: {found}")
    else:
        print("AutoDock-GPU compilation finished but no binary was found in bin/.")
        sys.exit(1)

    # make sure it’s executable (should already be, but belt & suspenders)
    st = os.stat(found)
    os.chmod(found, st.st_mode | stat.S_IXUSR)
    print(f"[BUILD] AutoDock-GPU ready: {found}")

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


def main():
    CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
    print("=" * 50)
    print("Welcome to the Ultidock Setup")
    print("=" * 50)

    mode = _normalize_mode(input("Select run mode (GPU = NVidia-CUDA/OpenCL, AMD-OpenCL or CPU or auto) [default: auto]: ") or "auto")
    if mode in ("opencl", "cuda", "cpu"):
        GPU_TYPE = mode.upper()          
    else:  # 'gpu' or 'auto'
        GPU_TYPE = detect_gpu()          
        
    # Ask for directory paths
    LIGANDS_DIR = ask_for_input("Enter the path for ligand files", os.path.join(CURRENT_DIR, "LIGANDS_DIR"))
    DOCKING_DIR = ask_for_input("Enter the path for docking files", os.path.join(CURRENT_DIR, "DOCKING_DIR"))
    ANALYSIS_DIR = ask_for_input("Enter the path for analysis files", os.path.join(CURRENT_DIR, "ANALYSIS_DIR"))
    VINA_DIR = ask_for_input("Enter the path for Autodock Vina", os.path.join(CURRENT_DIR, "VINA_DIR"))
    AUTODOCK_GPU_DIR = ask_for_input("Enter the path for Autodock-GPU files", os.path.join(CURRENT_DIR, "AUTODOCK_GPU_DIR"))
    MACRO_MOL_DIR = ask_for_input("Enter the path for macro molecule of your choice", os.path.join(CURRENT_DIR, "MACRO_MOL_DIR"))
    RESULTS_DIR = ask_for_input("Enter the path for results files", os.path.join(CURRENT_DIR, "RESULTS_DIR"))

    # Ask for the .wget file location
    wget_file_path = ask_for_input("Enter the path to the .wget file", os.path.join(CURRENT_DIR, "ligands.wget"))

    # Create required directories (including results)
    create_directory_if_needed(LIGANDS_DIR)
    create_directory_if_needed(DOCKING_DIR)
    create_directory_if_needed(ANALYSIS_DIR)
    create_directory_if_needed(VINA_DIR)
    create_directory_if_needed(AUTODOCK_GPU_DIR)
    create_directory_if_needed(MACRO_MOL_DIR)
    create_directory_if_needed(RESULTS_DIR)

    # Detect and compile AutoDock-GPU if needed
    detect_and_compile_autodock_gpu(AUTODOCK_GPU_DIR, GPU_TYPE, NUMWI)

    # Save configuration to config.py (only declaring paths; DB file is not created here)
    with open(os.path.join(ROOT_DIR,"docking", "config.py"), 'w') as config_file:
        config_file.write('# config.py\n')
        config_file.write('# Auto-generated config.py\n')
        config_file.write('import os\n\n')
        # Use the directory where config.py is located as the base directory.
        config_file.write('BASE_DIR = os.path.abspath(os.path.dirname(__file__))\n\n')
        config_file.write('LIGANDS_DIR = os.path.join(BASE_DIR, "LIGANDS_DIR")\n')
        config_file.write('DOCKING_DIR = os.path.join(BASE_DIR, "DOCKING_DIR")\n')
        config_file.write('ANALYSIS_DIR = os.path.join(BASE_DIR, "ANALYSIS_DIR")\n')
        config_file.write('VINA_DIR = os.path.join(BASE_DIR, "VINA_DIR")\n')
        config_file.write('AUTODOCK_GPU_DIR = os.path.join(BASE_DIR, "AUTODOCK_GPU_DIR")\n')
        config_file.write('MACRO_MOL_DIR = os.path.join(BASE_DIR, "MACRO_MOL_DIR")\n')
        config_file.write('RESULTS_DIR = os.path.join(BASE_DIR, "RESULTS_DIR")\n')
        config_file.write('GPU_TYPE = "' + GPU_TYPE + '"\n')
        config_file.write('DB_PATH = os.path.join(RESULTS_DIR, "ultidock_results.db")\n')
        config_file.write(f'NUMWI = "{NUMWI}"\n')

        print("Configuration saved to config.py and default directories are ensured!")

    # Download the ligands using the URLs from the .wget file
    download_ligands_from_file(wget_file_path, LIGANDS_DIR)


if __name__ == "__main__":
    main()
