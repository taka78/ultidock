# Contributing

Discuss a scientific behavior change with a reproducible example and an explicit
expected outcome. Distinguish formatting repairs from chemical changes and site
selection changes from docking-engine changes.

## Local setup

```bash
git clone --branch gmx-dev https://github.com/taka78/ultidock.git
cd ultidock
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest
```

Native integrations have additional prerequisites; see [Testing](testing.md).
For MD development, install `.[dev,md]` and the external toolchain, then use
[MD validation](md-validation.md).
Keep generated grids, database files, local configuration and example workspaces
out of source control. Preserve small source fixtures needed for reproducible tests.

## Changes to review

Describe the concrete failure and resulting behavior. Include relevant test results
and any chemistry or compatibility implications. When changing pose selection,
verify that exported scores refer to the saved coordinates. When changing receptor
preparation, check heavy-atom preservation and provenance. For CLI changes, update
help, examples and the configuration reference together.

Add runnable datasets under `examples/<name>/` with an `example-run.py` entry point,
input provenance and a readable explanation. Use the shared staging helper to keep
runs independent. Do not present comparison molecules as benchmark actives/decoys
without supporting labels.

Before proposing documentation changes, run the strict Sphinx build described in
[Testing](testing.md). Keep instructions consistent with executable help and avoid
publishing illustrative output as if it came from a docking experiment.
