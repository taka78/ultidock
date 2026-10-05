# Ultidock Release Checklist

Use this before tagging a citable source-checkout release.

- Run `ultidock doctor`. done
- Run `ultidock example run quickstart`. done
- Run CLI smoke checks for `known-site`, `cavity`, `blind`, and `report`. done
- Run at least one CaV-EMPS DUD-E cavity-recovery smoke benchmark. done
- Run the normalized site-prediction evaluator on a tiny fixture with an empty method output. done
- Confirm reports include command, config, version, timestamp, and git commit. done
- Confirm generated benchmark datasets, maps, and run results are absent from the Git source archive. done
- Run the Python tests and check `git diff --check`. done
- Update `CHANGELOG.md` and `CITATION.cff` for the exact release version and date.
- Verify the release source archive can be installed from a checkout with `pip install -e .`. done
- Enable the public GitHub repository in Zenodo before publishing a GitHub release.
- Commit the candidate, tag it, and publish the GitHub release. Zenodo then archives that release and assigns a DOI.
- Check the Zenodo record and add its DOI to repository citation metadata after it exists.
