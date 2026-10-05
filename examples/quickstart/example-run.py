#!/usr/bin/env python3
"""An interactive Next/Back/Exit teacher for a real Ultidock example."""

import argparse
import importlib.util
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

EXAMPLES = Path(__file__).resolve().parent.parent
ROOT = EXAMPLES.parent


def tool(name, relative):
    local = ROOT / relative
    return str(local) if local.is_file() and os.access(local, os.X_OK) else shutil.which(name)


def environment(mode, site_method="cav-emps"):
    lines = []
    missing = []
    for module in ("numpy", "scipy", "psutil", "pandas"):
        found = importlib.util.find_spec(module) is not None
        lines.append(f"{'OK' if found else 'MISSING'}  Python: {module}")
        if not found:
            missing.append(f"Python package {module}; install the project dependencies")
    if mode == "cpu":
        found = tool("vina", "docking/VINA_DIR/bin/vina")
        lines.append(f"{'OK' if found else 'MISSING'}  Vina: {found or 'not installed'}")
        if not found:
            missing.append("AutoDock Vina; on Ubuntu: sudo apt install autodock-vina")
    grid = tool("autogrid4", "docking/AUTODOCK_GPU_DIR/autogrid/autogrid4")
    lines.append(f"AutoGrid: {grid or 'setup will try to build the bundled source'}")
    babel = shutil.which("obabel")
    lines.append(f"Open Babel: {babel or 'not installed'}")
    if not babel:
        missing.append("Open Babel for missing-hydrogen recovery; on Ubuntu: sudo apt install openbabel")
    if mode == "gpu":
        lines.append("GPU runtime/device validation happens in setup; GPU mode does not silently fall back to CPU.")
    if site_method != "cav-emps":
        relative = ("external/fpocket/bin/fpocket" if site_method == "fpocket"
                    else "external/bin/prank")
        local = ROOT / relative
        found = local.is_file() and os.access(local, os.X_OK)
        lines.append(f"{'OK' if found else 'MISSING'}  {site_method}: {local}")
        if not found:
            missing.append(f"{site_method}; from the application workspace: "
                           f"bash scripts/install_pocket_tools.sh {site_method}")
        if site_method == "p2rank":
            lines.append("P2Rank also needs a compatible Java runtime; see SETUP.md.")
    return lines, missing


def pages(args):
    inputs = ("6CM4-edited.pdbqt; haloperidol, escitalopram, morphine"
              if args.example == "d2-antipsychotics" else "5i6x_edited.pdbqt; escitalopram")
    checks, missing = environment(args.mode, args.site_method)
    return [
        ("Welcome", "You will run a real receptor/ligand example, not a simulated result.\n"
         "First we check the tools, then inspect the inputs, choose sites, and run docking.\n"
         f"Example: {args.example}   Engine: {args.mode}   Sites: {args.site_method}\n"
         "Enter = Next, b = Back, q = Exit. Nothing runs until the final Run screen."),
        ("Check your installation", "\n".join(checks) + "\n\n" +
         ("Install the missing items, then restart this wizard:\n- " + "\n- ".join(missing)
          if missing else "Python dependencies and the checked tools are available.") +
         "\nFull native build instructions: SETUP.md. This teacher does not run sudo or install packages."),
        ("Understand the inputs", f"{inputs}\n\n"
         "The receptor is the target structure. Ligand PDBQT files contain candidate molecules,\n"
         "charges and rotatable bonds. MolGuard checks formatting and prepares receptors.\n"
         "Missing donor hydrogens can trigger recovery; read any chemistry warnings.\n"
         "Only the example's listed inputs are copied into a fresh timestamped workspace."),
        ("Understand the binding sites", f"Selected method: {args.site_method}\n\n"
         "A docking box defines where the engine searches. CaV-EMPS uses receptor geometry\n"
         "and AutoGrid maps; fpocket and P2Rank are alternative installed pocket finders.\n"
         "Site IDs such as S1 are identifiers, not proof that a site is correct.\n"
         "The run will prepare grids and dock each listed ligand at the proposed sites."),
        ("Run docking", "Next starts preparation, site discovery and docking. This may take time.\n"
         "The terminal will show the workspace and progress. Outputs stay in that workspace.\n"
         "After completion, open RESULTS_DIR's CSV and follow docking_file to the best\n"
         "PDBQT pose. A score ranks poses; it does not establish measured affinity or activity."),
    ], missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", choices=("d2-antipsychotics", "sert-escitalopram"), default="d2-antipsychotics")
    parser.add_argument("--mode", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--site-method", choices=("cav-emps", "fpocket", "p2rank"), default="cav-emps")
    parser.add_argument("--dry-run", action="store_true", help="Read all lesson screens without prompts or docking.")
    parser.add_argument("--yes", action="store_true", help="Explicitly run all stages without prompts (for scripts).")
    args = parser.parse_args()
    screens, missing = pages(args)
    if not args.dry_run and not args.yes and not sys.stdin.isatty():
        parser.error("Interactive teacher needs a terminal. Use --dry-run to read it or --yes to run explicitly.")
    index = 0
    while index < len(screens):
        title, body = screens[index]
        print(f"\n{'=' * 60}\nStep {index + 1} of {len(screens)}: {title}\n{'=' * 60}\n{body}\n", flush=True)
        if args.dry_run or args.yes:
            answer = ""
        else:
            try:
                answer = input("[Enter] Next   [b] Back   [q] Exit: ").strip().lower()
            except (EOFError, KeyboardInterrupt):
                print("\nExited. No docking was started.")
                return
        if answer in {"q", "quit", "exit"}:
            print("Exited. No docking was started.")
            return
        if answer in {"b", "back"}:
            index = max(0, index - 1)
            continue
        if answer not in {"", "n", "next"}:
            print("Use Enter, b, or q.")
            continue
        if index == 1 and missing and not args.dry_run:
            parser.exit(1, "Resolve the missing prerequisites above, then restart the wizard.\n")
        index += 1
    command = [sys.executable, "-u", str(EXAMPLES / args.example / "example-run.py"),
               args.site_method, "--mode", args.mode]
    print("Pipeline command: " + shlex.join(command), flush=True)
    if args.dry_run:
        print("Preview complete. No files were staged and no docking was started.")
        return
    result = subprocess.run(command, cwd=EXAMPLES / args.example)
    if result.returncode:
        parser.exit(result.returncode, "The example did not complete. Read the error and workspace logs above.\n")
    print("\nFinished. Follow the printed workspace path to RESULTS_DIR and DOCKING_DIR.\n"
          "Compare poses in a molecular viewer and retain preparation warnings with your results.\n"
          "Continue with docs/source/user-guide/results-reports.md for interpretation.")


if __name__ == "__main__":
    main()
