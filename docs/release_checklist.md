# Ultidock Release Checklist

Use this before tagging a citable source-checkout release.

- Run `ultidock doctor`.
- Run `ultidock example run quickstart`.
- Run CLI smoke checks for `known-site`, `cavity`, `blind`, and `report`.
- Run at least one CaV-EMPS DUD-E cavity-recovery smoke benchmark.
- Run the normalized site-prediction evaluator on a tiny fixture with an empty method output.
- Confirm reports include command, config, version, timestamp, and git commit.
- Confirm generated benchmark datasets, maps, and run results are absent from the Git source archive.
- Run the Python tests and check `git diff --check`.
- Update `CHANGELOG.md` and `CITATION.cff` for the exact release version and date.
- Verify the release source archive can be installed from a checkout with `pip install -e .`.
- Enable the public GitHub repository in Zenodo before publishing a GitHub release.
- Commit the candidate, tag it, and publish the GitHub release. Zenodo then archives that release and assigns a DOI.
- Check the Zenodo record and add its DOI to repository citation metadata after it exists.

This release uses the documented source-checkout installation path. Do not attach
or advertise the current Python wheel: it omits workflow directories and its
`ultidock` entry point fails outside a checkout. Fix and test wheel packaging
before distributing one. Keep large benchmark datasets and paper-run evidence
as a separate data archive if they need their own DOI.
