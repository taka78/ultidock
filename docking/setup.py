
import os
import subprocess
import sys
import shutil

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
    if shutil.which("nvidia-smi"):
        try:
            subprocess.run(["nvidia-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            return "NVIDIA"
        except subprocess.CalledProcessError:
            print("NVIDIA GPU detected but nvidia-smi failed to execute.")
            return "Unknown"

    # Check for AMD GPU using ROCm tools
    if shutil.which("rocm-smi"):
        try:
            subprocess.run(["rocm-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            return "AMD"
        except subprocess.CalledProcessError:
            print("AMD GPU detected but rocm-smi failed to execute.")
            return "Unknown"

    # Default to CPU if no GPU is detected
    print("No compatible GPU detected. Using CPU mode.")
    return "CPU"


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


def main():
    CURRENT_DIR = os.getcwd()

    print("Welcome to the Ultidock Setup")

    GPU_TYPE = detect_gpu()  # This won't crash the script anymore

    # Ask for directory paths
    LIGANDS_DIR = ask_for_input("Enter the path for ligand files", os.path.join(CURRENT_DIR, "LIGANDS_DIR"))
    DOCKING_DIR = ask_for_input("Enter the path for docking files", os.path.join(CURRENT_DIR, "DOCKING_DIR"))
    ANALYSIS_DIR = ask_for_input("Enter the path for analysis files", os.path.join(CURRENT_DIR, "ANALYSIS_DIR"))
    VINA_GPU_DIR = ask_for_input("Enter the path for Vina-GPU-2.0 files", os.path.join(CURRENT_DIR, "VINA_GPU_DIR"))
    VINA_DIR = ask_for_input("Enter the path for Autodock Vina", os.path.join(CURRENT_DIR, "VINA_DIR"))
    MACRO_MOL_DIR = ask_for_input("Enter the path for macro molecule of your choice.", os.path.join(CURRENT_DIR, "MACRO_MOL_DIR"))

    # Ask for the .wget file location
    wget_file_path = ask_for_input("Enter the path to the .wget file", os.path.join(CURRENT_DIR, "ligands.wget"))

    # Detect GPU and select the appropriate docking engine

    # Save these to config.py
    with open('config.py', 'w') as config_file:
        config_file.write(f'# config.py \n')
        config_file.write(f'# Paths for directories \n')
        config_file.write(f'LIGANDS_DIR = "{LIGANDS_DIR}"\n')
        config_file.write(f'DOCKING_DIR = "{DOCKING_DIR}"\n')
        config_file.write(f'ANALYSIS_DIR = "{ANALYSIS_DIR}"\n')
        config_file.write(f'VINA_DIR = "{VINA_DIR}"\n')
        config_file.write(f'VINA_GPU_DIR = "{VINA_GPU_DIR}"\n')  # Fixed overwrite issue
        config_file.write(f'MACRO_MOL_DIR = "{MACRO_MOL_DIR}"\n')
        config_file.write(f'GPU_TYPE = "{GPU_TYPE}"\n')  # Save GPU detection result

    # Create required directories
    create_directory_if_needed(LIGANDS_DIR)
    create_directory_if_needed(DOCKING_DIR)
    create_directory_if_needed(ANALYSIS_DIR)
    create_directory_if_needed(VINA_GPU_DIR)
    create_directory_if_needed(VINA_DIR)
    create_directory_if_needed(MACRO_MOL_DIR)

    print("Configuration saved to config.py and default directory files created!")

    # Download the ligands using the URLs from the .wget file
    download_ligands_from_file(wget_file_path, LIGANDS_DIR)


if __name__ == "__main__":
    main()
