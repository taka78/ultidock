<h1 align=\"center\">Ultidock: High-Throughput Docking Pipeline</h1>

<p align=\"center\">
  <div style="text-align: center;">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/taka78/ultidock/dev-beta/traffic-badge.json" alt="GitHub Traffic Badge" />
</p>

<h2>📖 Overview</h2>
<p>
Ultidock is a powerful and fully automated molecular docking pipeline built around AutoDock Vina, optimized for large-scale ligand screening with SQLite database integration for efficient data handling.
</p>

<h2>✨ Features</h2>
<ul>
  <li><strong>Fully Automated Workflow:</strong> Setup, extraction, docking, and analysis.</li>
  <li><strong>Optimized Performance:</strong> Multithreaded docking processes.</li>
  <li><strong>SQLite Database Integration:</strong> Structured data storage for large datasets.</li>
  <li><strong>Flexible Analysis:</strong> Advanced SQL-based data filtering and export.</li>
  <li><strong>User-Friendly Config:</strong> Simple configuration through a single <code>config.py</code>.</li>
</ul>

<h2>⚙️ Requirements</h2>
<ul>
  <li>Python 3.10+</li>
  <li>AutoDock Vina (included)</li>
  <li>SQLite3 (built-in Python)</li>
  <li>Pandas, NumPy</li>
</ul>

<h2>🚀 Quick Start</h2>
<h3>1. Clone the repository:</h3>

<pre><code>git clone https://github.com/taka78/ultidock.git
cd ultidock
</code></pre>

<h3>2. Adjust the expected analysing standarts for your molecule from <code>analyse_docking_results.py</code>:</h3>
<pre><code>
    DEFAULT_AFFINITY = -7.0        # kcal/mol // you should change this according to how much chemically active your macromolecule.
    DEFAULT_RMSD_LB = 5.0          # Å // you should change this according to how big your macromolecule's docking site is.
    DEFAULT_RMSD_UB = 10.0         # Å // you should change this according to how big your macromolecule's docking site is.
    DEFAULT_MIN_MODEL = 2          # integer // you should change this according to how picky you are.
</code></pre>


<h3>3. Run the pipeline:</h3>
<pre><code>python3 /path/to/your/ultidock/docking/run.py
</code></pre>


<h2>🛠 Configuration (<code>config.py</code>)</h2>
<p>Defined all paths and default thresholds in a single configuration file for seamless adjustments automatically when you run run.py:</p>

<pre><code>BASE_DIR = '/path/to/your/ultidock'
DB_PATH = f"{BASE_DIR}/results/ultidock_results.db"
LIGANDS_DIR = f"{BASE_DIR}/docking/LIGANDS_DIR"
MACRO_MOL_DIR = f"{BASE_DIR}/docking/MACRO_MOL_DIR"
DOCKING_DIR = f"{BASE_DIR}/docking"
VINA_DIR = f"{BASE_DIR}/docking"
</code></pre>
<p>Edit <code>config.py</code> to set paths and parameters according to your needs.</p>

<h2>🔍 Results</h2>
<p>Docking results are stored in SQLite (<code>results/ultidock_results.db</code>). Results will be exported to CSV automatically.</p>

<h2>💾 Exporting Results</h2>

<pre><code>
# Export to Excel (if you need an .xlsx file for whatever the reason.)
python docking/analyse_docking_results.py --out results.xlsx
</code></pre>

<h2>🤝 Contributing</h2>
<p>Contributions are welcome! Please open an issue or submit a pull request.</p>
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


  <hr />

<p align=\"center\">⭐ If you find Ultidock useful, please star the repository!</p>

</body>
</html>
