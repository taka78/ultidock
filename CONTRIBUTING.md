# Contributing to Ultidock

Bug reports, support questions, documentation fixes, and code contributions are welcome.

## Issues and support

Search [existing issues](https://github.com/taka78/ultidock/issues) before opening one.
For bugs, include your Ultidock version, OS, Python version, command, and traceback.
Describe expected versus actual behavior and attach a small reproducible input.
For usage help, open an issue describing your goal and what you have tried.
Discuss changes to scientific behavior in an issue before implementing them.

## Development and tests

Fork the repository, clone your fork, and create a branch for your change.
From the repository root, use Python 3.10+ and install the development dependencies:

```bash
python -m pip install --upgrade pip 'setuptools>=68' wheel
python -m pip install -e '.[dev]'
python -m pytest molguard/tests/ test_dude_benchmark_helpers.py
```

CPU tests run in GitHub Actions without a GPU; see [Testing](docs/source/development/testing.md) for more checks.

## Pull requests

Submit a PR against the branch you started from; link the issue and explain the change.
Include relevant regression tests, updated documentation, and the test results.
Keep generated grids, databases, local configuration, and run outputs out of the PR.
