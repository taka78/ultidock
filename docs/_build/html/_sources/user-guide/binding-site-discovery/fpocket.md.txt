# fpocket

fpocket is an external geometric pocket finder. Ultidock's adapter runs the local
tool, reads its ranked pocket output and generates a centers TSV for docking.
Install it from the checkout, then inspect tool discovery:

```bash
bash scripts/install_pocket_tools.sh fpocket
ultidock doctor
```

Run the D2 example with fpocket sites:

```bash
ultidock example run d2-antipsychotics fpocket --mode cpu
```

Generate boxes without launching docking:

```bash
ultidock pocket-box --method fpocket --receptor /absolute/receptor.pdbqt \
  --output /absolute/sites.tsv --autosites 6 --box-size 35
```

Use `--tool /path/to/fpocket` to select a specific executable. Preserve that
version and its raw pocket outputs when comparing against CaV-EMPS or P2Rank.
The wrapper's box dimensions and maximum site count are part of the experiment,
not just presentation settings. Inspect boxes for alignment with the receptor's
coordinate frame before docking.

For a combined custom workflow, `ultidock fpocket --help` shows the receptor,
box and output options; ordinary pipeline options are forwarded. Its `--dry-run`
can still execute pocket prediction to create boxes—it is not a no-work preview.
See [CLI reference](../../reference/cli.md).
