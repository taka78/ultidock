# Interactive first-run teacher

```bash
ultidock example run quickstart
```

This opens five short screens: welcome, installation check, receptor/ligand
preparation, binding sites, and running/reading docking results. Press **Enter**
to go Next, **b** to go Back, or **q** to exit. Docking starts only when you advance
past the final Run screen. Missing requirements are explained with installation
hints; the teacher does not install packages or invoke sudo.

The default lesson runs the real D2 example with haloperidol, escitalopram and
morphine, using CPU Vina and CaV-EMPS. Choose another supported lesson or engine:

```bash
ultidock example run quickstart --example sert-escitalopram
ultidock example run quickstart --mode gpu
ultidock example run quickstart --site-method fpocket
```

Read every screen without running anything:

```bash
ultidock example run quickstart --dry-run
```

For an explicitly unattended run, use `--yes`. Without it, non-interactive input
is rejected rather than silently starting a long calculation. The selected example
prints its own timestamped workspace; old outputs are preserved. This replaces the
old artifact generator and its synthetic scores. No fixture report is generated.

Start with [Installation](../../docs/source/getting-started/installation.md), then
[First Docking Run](../../docs/source/getting-started/first-docking-run.md). Source-checkout
users can replace `ultidock` with `/usr/bin/python3 -m cli.ultidock` from the repo root.
