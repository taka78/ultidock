# Quick Start

The quickstart is an interactive teacher, not an artificial result generator.
Launch it in a terminal:

```bash
ultidock example run quickstart
```

At each screen, press **Enter** for Next, **b** for Back, or **q** to exit.

1. **Welcome:** see the example, engine and site method.
2. **Installation:** check Python dependencies and native tools. Missing tools
   produce installation hints; the teacher does not run privileged installers.
3. **Inputs:** learn what the receptor and ligand PDBQT files represent.
4. **Sites:** understand the search boxes and the selected pocket method.
5. **Run:** review the calculation and press Next to launch the real pipeline.

Nothing is staged and no docking starts until you advance past the final Run
screen. The default lesson uses D2/6CM4 with haloperidol, escitalopram and morphine,
CPU Vina, and CaV-EMPS. The example then prints a new timestamped workspace.

Read the complete lesson without prompts or execution:

```bash
ultidock example run quickstart --dry-run
```

Choose a different lesson or backend at launch:

```bash
ultidock example run quickstart --example sert-escitalopram
ultidock example run quickstart --mode gpu
ultidock example run quickstart --site-method fpocket
```

Use `--yes` only when you intend to run unattended. Without that flag, a
non-interactive terminal is rejected instead of silently starting docking.
For a direct command and the expected files, continue to [First Docking Run](first-docking-run.md).
