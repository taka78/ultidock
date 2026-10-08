# Unknown binding-site workflow

Use automatic site discovery when you do not have a defensible manual binding box.
Put a receptor and ligand library in the default workspace folders described in
[Dock your own molecules](../user-guide/start-docking.md). The pipeline prepares
the receptor and selects the available GPU or CPU backend automatically.

## 1. Choose the inputs

Place receptor `.pdb`, `.mol2` or `.pdbqt` files in `docking/MACRO_MOL_DIR/`.
Place prepared ligand PDBQTs in `docking/LIGANDS_DIR/`, or use
`docking/ligands.wget` to download archives. Keep only the molecules intended
for this run. Do not leave a reference ligand embedded in a receptor used for
receptor-only site prediction.

## 2. Start docking

```bash
ultidock cavity
```

This prepares receptors, proposes CaV-EMPS sites and docks ligands. If you
placed local ligands in `LIGANDS_DIR`, add `--skip-wget` to exclude the bundled
example download. To preview the generated command and run configuration
without docking, add `--dry-run`.

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
