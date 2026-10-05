# Testing

Install development dependencies, then run from the repository root:

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

Tests live under `molguard/tests/`, including coverage for the broader workflow,
CLI and packaging. Some integrations depend on installed native tools; review skip
reasons and failures rather than treating a partial environment as a complete
end-to-end validation.

## Focused checks

```bash
python -m pytest molguard/tests/test_receptor_prep.py
python -m pytest molguard/tests/test_example_site_methods.py molguard/tests/test_quickstart.py
python -m pytest molguard/tests/test_package_installation.py
python -m cli.ultidock example run quickstart --dry-run
python -m cli.ultidock example run d2-antipsychotics --dry-run
```

The packaging test builds and installs artifacts in a temporary environment and
runs commands outside the checkout. It checks that runtime resources are present
and generated files stay in the managed workspace.

Dry runs and mocked dispatch tests verify control flow, not docking accuracy.
For native smoke validation, run the D2 example with installed tools, inspect the
workspace and verify that every exported score is associated with a real output
pose. Benchmark scientific behavior separately on declared reference datasets.

## Documentation build

```bash
python -m pip install -r docs/requirements.txt
make -C docs html SPHINXOPTS="-n -W --keep-going"
```

Open `docs/build/html/index.html`. The strict build treats broken internal
references and other Sphinx warnings as failures. Read the Docs uses the same
configuration and pinned documentation dependencies, without building native engines.

The tutorial-template layout keeps `conf.py`, `index.rst` and guide pages in
`docs/source/`, with generated output in `docs/build/`. Without Make, the
equivalent command is:

```bash
python -m sphinx -b html -n -W --keep-going docs/source docs/build/html
```
