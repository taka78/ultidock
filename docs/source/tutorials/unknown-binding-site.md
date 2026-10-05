# Unknown binding-site workflow

Use automatic site discovery when you do not have a defensible manual binding box.
Prepare one receptor and a small ligand set before scaling up. The example paths
below are placeholders for your own absolute directories.

## 1. Prepare a dedicated input directory

Keep only the intended receptor in `/absolute/project/receptors` and prepared ligand
PDBQTs in `/absolute/project/ligands`. Review [receptor preparation](../user-guide/preparing-receptor.md)
and [ligand preparation](../user-guide/preparing-ligands.md). Do not leave a reference
ligand embedded in a receptor used for receptor-only site prediction.

## 2. Preview the run

```bash
ultidock cavity --autosites 6 --mode cpu --skip-wget \
  --macro-mol-dir /absolute/project/receptors \
  --ligands-dir /absolute/project/ligands \
  --output-dir /absolute/project/run-cavity --dry-run
```

This writes run configuration and prints the pipeline command; it does not validate
chemical correctness or prove that native tools can execute. Inspect the paths,
then repeat without `--dry-run` to run preparation, site discovery and docking.

## 3. Inspect site coverage

Load the prepared receptor, generated centers/boxes and representative best poses.
A high-scoring pose outside your biologically relevant region may indicate a site
selection problem. Increase the site budget only with a clear evaluation plan;
more sites cost more docking jobs.

## 4. Compare an alternative method

Use [fpocket](../user-guide/binding-site-discovery/fpocket.md) or
[P2Rank](../user-guide/binding-site-discovery/p2rank.md) with the same receptor and
box budget. Keep their prediction artifacts and separate output directories. If
experimental reference sites become available, evaluate localization independently
from pose scores, as described in [Binding-site prediction](../scientific-background/binding-site-prediction.md).
