<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Ultidock – Beta Branch</title>
</head>
<body style="font-family: Arial, sans-serif; line-height: 1.6; margin: 2rem; background-color: #fdfdfd; color: #333;">

  <h1>Ultidock Project – Beta Channel</h1>
  <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/taka78/ultidock/dev-beta/traffic-badge.json" alt="GitHub Traffic Badge" />

  <p>Ultidock is a parallelized, semi-automated docking pipeline designed for high-throughput molecular screening using AutoDock Vina. This <strong>beta branch</strong> is stable and functional, built to simplify ligand screening for researchers and students with moderate hardware setups.</p>

  <hr />

  <h2>How It Works</h2>
  <p>This version of Ultidock performs the following tasks:</p>
  <ul>
    <li>Downloads ligand structures using a <code>.wget</code> list</li>
    <li>Prepares a grid around a supplied macromolecule</li>
    <li>Runs docking jobs in batches via AutoDock Vina</li>
    <li>Outputs structured docking results into a specified directory</li>
  </ul>

  <hr />

  <h2>Step-by-Step Instructions</h2>

  <h3>1. Prepare Your Macromolecule</h3>
  <ul>
    <li>Download your macromolecule from <a href="https://www.rcsb.org/">RCSB PDB</a> or <a href="https://alphafold.ebi.ac.uk/">AlphaFold</a>.</li>
    <li>Remove all water molecules and unnecessary chains.</li>
    <li><strong>Keep REMARK headers</strong> as shown in the example files.</li>
    <li>Convert it to <code>.pdbqt</code> format and place it in: <code>docking/MACRO_MOL_DIR/</code></li>
  </ul>

  <h3>2. Prepare the Ligand Download List</h3>
  <ul>
    <li>Create a file named <code>ligands.wget</code> that contains <code>wget</code> commands for your desired ligands.</li>
    <li>Example sources: <a href="https://zinc15.docking.org/">ZINC15</a>, <a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a>, <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>, <a href="https://go.drugbank.com/">DrugBank</a></li>
  </ul>

  <pre><code>wget http://some-ligand-url.com/ligand1.mol2
wget http://some-ligand-url.com/ligand2.mol2</code></pre>

  <p>Save this file in your main working directory: <code>docking/ligands.wget</code></p>

  <h3>3. Run the Docking Pipeline</h3>
  <p>From within the <code>docking/</code> directory, run:</p>

  <pre><code>python run.py</code></pre>

  <p>This will:</p>
  <ul>
    <li>Download ligands</li>
    <li>Generate a centered grid box</li>
    <li>Run docking jobs sequentially using multiple threads</li>
    <li>Save the output <code>.pdbqt</code> files in <code>DOCKING_DIR</code></li>
  </ul>

  <h3>4. Analyze Results</h3>
  <p>Once docking is complete, you can sort and analyze the output:</p>

  <pre><code>python output-analyses.py</code></pre>

  <p>This script parses <code>.pdbqt</code> results and outputs CSVs using Pandas. SQL-based analysis is planned in future versions.</p>

  <hr />

  <h2>Performance Notes</h2>
  <p>This version is designed for CPU-based processing. An NVMe SSD is recommended for best performance, but Optane or other high-end storage is no longer required thanks to simplified I/O handling.</p>

  <blockquote><strong>Warning:</strong> This tool is CPU-intensive. Ensure your system is cooled and monitored during batch runs.</blockquote>

  <hr />

  <h2>System Configuration (Used in Development)</h2>
  <ul>
    <li>CPU: AMD Ryzen 5 3600X</li>
    <li>RAM: 24 GB DDR4</li>
    <li>Storage: 1 TB NVMe SSD</li>
    <li>Environment: WSL + Python 3.10</li>
  </ul>

  <p>Running ~1.2 million ligands against a macromolecule (4H10) took about 3 days and generated ~80 GB of data.</p>

  <hr />

  <h2>What's Coming in dev-beta</h2>

  <p>The <code>dev-beta</code> branch improves on this version with:</p>
  <ul>
    <li>SQLite-based result logging</li>
    <li>Multi-thread-safe database inserts</li>
    <li>GPU acceleration (coming soon)</li>
    <li>Automated grid generation via <code>setup.py</code></li>
    <li>Cleaner, portable directory structure</li>
  </ul>

  <p>If you're looking for better performance tracking and large-scale scalability, try the <code>dev-beta</code> branch instead.</p>

  <hr />

  <h3>Disclaimer</h3>
  <p>This is a beta release. Features are evolving. Stability is not guaranteed. Use at your own risk.</p>

</body>
</html>
