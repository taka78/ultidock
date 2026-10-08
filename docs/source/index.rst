Ultidock Documentation
======================

Ultidock prepares molecular inputs, proposes binding sites, runs docking engines,
and keeps scores connected to their output poses. Follow :doc:`user-guide/start-docking`
to dock your own molecules, or use :doc:`tutorials/examples` to learn from
bundled inputs.

These pages describe the **gmx-dev branch**, including its optional GROMACS
molecular-dynamics workflow. For docking through simulation, start with
:doc:`getting-started/md-installation` and
:doc:`user-guide/molecular-dynamics/index`.

**New here?** Check :doc:`getting-started/requirements`, follow
:doc:`getting-started/installation`, then use :doc:`user-guide/start-docking`.

For a lightweight tour of the reports, use :doc:`getting-started/quick-start`.
For a native GPU docking example, follow :doc:`getting-started/first-docking-run`.
For a selected ZINC library or your own prepared ligands, use
:doc:`tutorials/virtual-screening`.
For many receptors against a ligand library, use :doc:`user-guide/hpc-batch-screening`.
For published method comparisons, use :doc:`benchmarks/index`.

.. toctree::
   :caption: Getting Started
   :maxdepth: 1

   getting-started/what-is-ultidock
   getting-started/requirements
   getting-started/installation
   getting-started/md-installation
   getting-started/quick-start
   getting-started/first-docking-run

.. toctree::
   :caption: User Guide
   :maxdepth: 2

   user-guide/start-docking
   user-guide/preparing-receptor
   user-guide/preparing-ligands
   user-guide/binding-site-discovery/index
   user-guide/docking-engines/index
   tutorials/examples
   tutorials/small-molecule-docking
   tutorials/unknown-binding-site
   tutorials/virtual-screening
   user-guide/hpc-batch-screening
   user-guide/results-reports
   user-guide/molecular-dynamics/index

.. toctree::
   :caption: Scientific Background
   :maxdepth: 1

   scientific-background/docking-workflow
   scientific-background/binding-site-prediction
   scientific-background/cav-emps-methodology
   scientific-background/scoring-ranking
   scientific-background/molecular-dynamics
   scientific-background/limitations

.. toctree::
   :caption: Benchmarks
   :maxdepth: 1

   benchmarks/index

.. toctree::
   :caption: Reference
   :maxdepth: 1

   reference/cli
   reference/configuration
   reference/file-formats
   reference/troubleshooting
   reference/md-cli
   reference/md-protocol
   reference/md-troubleshooting

.. toctree::
   :caption: Development
   :maxdepth: 1

   development/architecture
   development/contributing
   development/testing
   development/md-validation
   development/citation
