<h1 align=\"center\">Ultidock: High-Throughput Docking Pipeline</h1>

<p align=\"center\">
  <div style="text-align: center;">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/taka78/ultidock/dev-beta/traffic-badge.json" alt="GitHub Traffic Badge" />
</p>

<h2>📖 Overview</h2>
<p>
Ultidock is a fully automated and lightweight molecular docking pipeline built around AutoDock Vina. Designed for large-scale virtual screening, it streamlines every stage—from ligand preparation to result analysis—while maintaining flexibility and performance. With integrated SQLite support and a multithreaded architecture, Ultidock is ideal for users aiming to dock thousands of ligands without hassle.
</p>

<h2>✨ Core Features</h2>
<ul>
  <li><strong>End-to-End Automation:</strong> Ligand preparation, docking, scoring, and result analysis in one pipeline.</li>
  <li><strong>Multithreaded Performance:</strong> Utilizes all available CPU cores for fast, parallel docking.</li>
  <li><strong>Integrated SQLite Database:</strong> Efficient and structured result storage, ideal for high-throughput workflows.</li>
  <li><strong>Configurable Analysis:</strong> Easily filter and export data based on RMSD, binding energy, or your own criteria.</li>
  <li><strong>Minimal Setup:</strong> Tweak your workflow with a single <code>config.py</code>—no GUI required.</li>
</ul>

<h2>⚙️ Requirements</h2>
<p>Ultidock is designed to be simple to deploy with only essential dependencies:</p>
<ul>
  <li>Python 3.10 or newer</li>
  <li><a href="https://vina.scripps.edu/">AutoDock Vina</a> (included or preconfigured)</li>
  <li>SQLite3 (included with Python)</li>
  <li>Pandas, NumPy (install via <code>pip install -r requirements.txt</code>)</li>
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
<p>
Ultidock uses a central configuration file to define all relevant paths and default parameters. When you launch the pipeline, these values are loaded automatically, allowing you to focus on your molecules instead of managing folders.
</p>

<pre><code>BASE_DIR = '/your/path/to/ultidock'
DB_PATH = f"{BASE_DIR}/results/ultidock_results.db"
LIGANDS_DIR = f"{BASE_DIR}/docking/LIGANDS_DIR"
MACRO_MOL_DIR = f"{BASE_DIR}/docking/MACRO_MOL_DIR"
DOCKING_DIR = f"{BASE_DIR}/docking"
VINA_DIR = f"{BASE_DIR}/docking"
</code></pre>

<p>
Macromolecule zoning and docking regions are predicted automatically by analyzing the geometry of each ligand. The system calculates grid centers and sizes based on ligand spatial distribution, then generates all required <code>.gpf</code> and <code>.glg</code> files for <code>autogrid4</code> behind the scenes—no need to manually define binding sites or run external preparation tools.
</p>

---

<h2>🔍 Results</h2>
<p>
All docking outcomes are logged directly into an SQLite database located at <code>results/ultidock_results.db</code>. This structured format allows for fast queries, filtering, and large-scale result aggregation.
</p>

<p>
In addition to database storage, successful docking results are tracked by ligand filename. This makes it easy to manually inspect or reprocess specific ligand–macromolecule pairs. If a file encounters an error or is skipped, it won’t be silently dropped—you’ll know.
</p>

<p>
For quick access or spreadsheet compatibility, docking summaries are also exported to <code>.csv</code> files alongside the database. This gives you flexibility to review, visualize, or integrate results into your own analysis workflows without needing a database viewer.
</p>

---

<h2>💾 Exporting Results</h2>
<p>
To export results to Excel or other formats, you can use the built-in analysis script:
</p>

<pre><code># Export to Excel (.xlsx)
python docking/analyse_docking_results.py --out results.xlsx
</code></pre>

<p>
You can also modify this script to change filter logic (e.g., affinity thresholds, pose count), or to output in CSV, TSV, or other formats. All exports are based on the data already stored in the SQLite database for reliability.
</p>

<p>
You can customize output columns, filter criteria, or data formats with a few edits to the analysis script.
</p>

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

<h2>📎 Acknowledgements</h2>
  <p>
  Ultidock relies on the robust and widely used <a href="http://vina.scripps.edu">AutoDock Vina</a> software for molecular docking. If you use Ultidock, please also cite the original Vina publication:
  </p>

  <blockquote>
  Trott, O., & Olson, A. J. (2010). <em>AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.</em> Journal of Computational Chemistry, 31(2), 455–461.  
  <a href="https://doi.org/10.1002/jcc.21334">https://doi.org/10.1002/jcc.21334</a>
  </blockquote>
  
  <hr />

<p align=\"center\">⭐ If you find Ultidock useful, please star the repository!</p>

</body>
</html>
