import os
import subprocess
import sys
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

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
            return "NVIDIA"
        except subprocess.CalledProcessError:
            print("NVIDIA GPU detected but inaccessible.")

    if has_amd:
        try:
            subprocess.run(["rocm-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            print("AMD GPU detected.")
            return "AMD"
        except subprocess.CalledProcessError:
            print("AMD GPU detected but inaccessible.")

    print("No compatible GPU detected. Using CPU mode.")
    return "CPU"

def _normalize_mode(s: str) -> str:
    s = (s or "").strip().lower()
    return s if s in ("auto", "gpu", "cpu") else "auto"


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

def detect_and_compile_autodock_gpu(AUTODOCK_GPU_DIR, GPU_TYPE):
    if GPU_TYPE in ("NVIDIA", "AMD"):
        # Ensure AutoDock-GPU is compiled if needed
        # detect_and_compile_autodock_gpu(AUTODOCK_GPU_DIR, GPU_TYPE)
        compiler_script = os.path.join(SCRIPT_DIR, "autodock-gpu-compiler.sh")
        print(f"AutoDock-GPU is set up for {GPU_TYPE} GPU.")
        try:
            subprocess.run(["bash", compiler_script, AUTODOCK_GPU_DIR], check=True, cwd=AUTODOCK_GPU_DIR)
        except subprocess.CalledProcessError as e:
            print("Compiler script failed.")
            sys.exit(1)
    else:
        print("AutoDock will run in CPU mode.")

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

    mode = _normalize_mode(input("Select run mode (auto/gpu/cpu) [default: auto]: ") or "auto")
    if mode == "cpu":
        print("Forcing CPU mode.")
        GPU_TYPE = "CPU"
    elif mode == "gpu":
        # Try to detect which GPU vendor is present; if none, warn and fall back to CPU
        detected = detect_gpu()  # returns "NVIDIA", "AMD", or "CPU"
        if detected in ("NVIDIA", "AMD"):
            GPU_TYPE = detected
        else:
            print("Warning: No compatible GPU detected; continuing in CPU mode.")
            GPU_TYPE = "CPU"
    else:
        # auto
        GPU_TYPE = detect_gpu()  # returns "NVIDIA", "AMD", or "CPU"
        
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
    detect_and_compile_autodock_gpu(AUTODOCK_GPU_DIR, GPU_TYPE)

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

        print("Configuration saved to config.py and default directories are ensured!")

    # Download the ligands using the URLs from the .wget file
    download_ligands_from_file(wget_file_path, LIGANDS_DIR)


if __name__ == "__main__":
    main()
