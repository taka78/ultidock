<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
</head>
<body style="font-family: Arial, sans-serif; line-height: 1.6; margin: 2rem; background-color: #f8f8f8; color: #333;">

  <h1>Ultidock Project – Dev Channel</h1>

  <p>Ultidock is a parallelized, automated docking pipeline designed for high-throughput molecular screening. If you're familiar with docking workflows, this tool is designed to minimize manual steps and accelerate large-scale simulations.</p>

  <hr />

  <h2>What’s New in the Beta Channel</h2>

  After a long (and slightly painful) detour through text-based read/write workflows, this release brings a full migration to SQL. Logging is now handled through an SQLite database with multi-thread-safe access, drastically improving performance and scalability.</p>

  <ul>
    <li>Switched from flat file parsing to structured <strong>SQLite logging</strong></li>
    <li>Introduced <strong>multi-thread-safe database access</strong> with locking</li>
    <li>Database includes: ligand name, receptor, affinity, RMSD, and timestamp</li>
    <li>Setup is fully automated via <code>setup.py</code>, which generates <code>config.py</code> with clean paths</li>
    <li>Paths are now relative to the working directory for better portability</li>
    <li>Results are stored in a <code>results/</code> directory, with the database created on first run</li>
    <li>Foundation laid for future features: resumable runs, GPU acceleration, and dashboard integration</li>
  </ul>

  <hr />

  <h2>How to Use</h2>

  <ol>
    <li>Create a <code>wget</code> file containing ligand download commands.</li>
    <li>Place the file in your <code>docking/</code> directory.</li>
    <li>Add your macromolecule structure to the <code>MACRO_MOL_DIR</code>.</li>
    <li>Run the main script:<br /><code>python run.py</code></li>
    <li>Analyse the database for good candidates.<li>

  </ol>

  <p>This will:</p>
  <ul>
    <li>Download and preprocess ligand files</li>
    <li>Automatically calculate and center the docking grid</li>
    <li>Start the docking process using AutoDock Vina with threading</li>
    <li>Log results to the database during execution</li>
  </ul>

  <p>To analyze results, use:</p>
  <pre><code>python output-analyses.py</code></pre>

  <p>A direct SQL integration with the analysis tool is planned, but for now, this script uses Pandas to extract and export results, typically to CSV.</p>

  <hr />

  <h2>Performance Notes</h2>

  <p>Early versions of Ultidock were heavily reliant on sequential text I/O, which significantly slowed down large-scale runs. SQL support changes everything. Once your data reaches millions of ligands or gigabytes of processed results, having indexed and structured access becomes critical.</p>

  <p><strong>SQLite now handles all result logging.</strong> If you're using a standard NVMe SSD, performance should be excellent. There's no longer a need to consider specialty drives like Optane for I/O bottlenecks.</p>

  <blockquote><strong>Note:</strong> This script is computationally intensive. Ensure proper CPU cooling and monitor system resources during execution.</blockquote>

  <hr />

  <h2>System Configuration (Dev Setup)</h2>

  <ul>
    <li>CPU: AMD Ryzen 5 3600X</li>
    <li>RAM: 24 GB DDR4</li>
    <li>Storage: 1 TB NVMe SSD</li>
    <li>Environment: WSL (Linux)</li>
  </ul>

  <p>Running ~1.2 million ligands against a single macromolecule (4H10) took around 3 days, generating approximately 80 GB of data.</p>

  <hr />

  <h2>Collaboration & Future Work</h2>

  <p>If you have access to more powerful compute infrastructure or are interested in contributing, particularly toward GPU integration or AI-assisted result filtering, feel free to reach out.</p>

  <p>This project is being developed by a physicist—not a molecular biologist—but the goal is to make advanced simulations more accessible, modular, and fast. Future directions include dynamic simulations and adaptive docking methods.</p>

  <hr />

  <hr />

  <h2>📖 Citation</h2>
  
  <p>If you use <strong>Ultidock</strong> in your research, publication, or automated pipeline, please consider citing it as:</p>
  
  <blockquote>
    Turgut, T. (2025). <em>Ultidock: A Lightweight Parallelized Docking Pipeline for Ligand Screening</em>. GitHub Repository. 
    <a href="https://github.com/taka78/ultidock">https://github.com/taka78/ultidock</a>
  </blockquote>
  
  <p>
    You are free to use and modify this software under the MIT License. 
    However, citation and credit are appreciated to support continued development.
  </p>


  <h3>Disclaimer</h3>
  <p>This is a beta release. Features are evolving. Stability is not guaranteed. Use at your own risk.</p>

<div style="text-align: center;">
  <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/taka78/ultidock/dev-beta/traffic-badge.json" alt="GitHub Traffic Badge" />

</body>
</html>
