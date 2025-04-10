<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Ultidock Project – Beta</title>
</head>
<body style="font-family: Arial, sans-serif; line-height: 1.6; margin: 2rem; background-color: #f8f8f8; color: #333;">

  <h1>🧪 Ultidock Project – Beta Channel</h1>

  <p>If you're familiar with docking workflows, this project can be a powerful tool in your research. Everything is now automated and streamlined for multi-threaded high-throughput simulations.</p>

  <hr />

  <h2>🚀 What’s New in the Beta Channel?</h2>

  <p><strong>WE NO LONGER NEED INTEL OPTANE!</strong><br />
  After an embarrassing journey through file-based hell, the dev finally discovered why the rest of the world uses SQL.<br />
  Ultidock now logs docking results directly to an SQLite database, with safe concurrent access and blazing-fast reads.</p>

  <ul>
    <li>🔄 Switched from text parsing to structured <strong>SQLite-based logging</strong></li>
    <li>🧵 Safe <strong>multi-threaded database access</strong> with connection locking</li>
    <li>📦 DB logs: ligand name, receptor ID, affinity, RMSD, timestamp</li>
    <li>⚙️ Fully automated <code>setup.py</code> generates paths in <code>config.py</code></li>
    <li>📁 Portable file structure — no more weird working directory issues</li>
    <li>💾 Results saved to <code>results/</code> and DB is created automatically</li>
    <li>🔥 Ready for future: resumable runs, GPU acceleration, even web dashboards</li>
  </ul>

  <hr />

  <h2>📦 How to Use</h2>

  <ol>
    <li>Create a <code>wget</code> file with ligand download links.</li>
    <li>Place the file in your main <code>docking/</code> directory.</li>
    <li>Put your macromolecule file into the <code>MACRO_MOL_DIR</code>.</li>
    <li>Run the project with:<br />
      <code>python run.py</code>
    </li>
  </ol>

  <p>This will:</p>
  <ul>
    <li>Download and split ligand structures</li>
    <li>Auto-center grid for docking</li>
    <li>Launch Vina jobs in parallel</li>
    <li>Log output to a structured SQLite database</li>
  </ul>

  <p>Once docking is complete, you can analyze the results with:</p>
  <pre><code>python output-analyses.py</code></pre>

  <p>This script uses Pandas to filter and export results, optionally to CSV.</p>

  <hr />

  <h2>⚙️ Performance Notes</h2>

  <p>Previously, this tool struggled with disk I/O. It was slow. I blamed hardware. I even recommended Optane. I was wrong.</p>

  <p>Thanks to SQL integration, <strong>an NVMe SSD is more than enough now</strong>. SQLite handles reads/writes like a champ. But you still need a solid CPU and thermal solution.</p>

  <blockquote><strong>WARNING:</strong> This script is a legitimate CPU burner. If you’re running it on air cooling, you might wanna watch those temps.</blockquote>

  <hr />

  <h2>💻 My Setup</h2>

  <ul>
    <li>CPU: Ryzen 5 3600X</li>
    <li>RAM: 24 GB DDR4</li>
    <li>Storage: 1 TB NVMe SSD</li>
    <li>OS: Linux</li>
  </ul>

  <p>Docking 1.2 million ligands with a single macromolecule (4H10) took ~3 days and generated 80 GB of data.</p>

  <hr />

  <h2>🤝 Collaboration</h2>
  <p>If you’ve got powerful compute resources or want to help expand this project (e.g., adding GPU backend or AI-assisted filtering), I’d love to talk. I’m a physicist, not a bioinformatician — I’m just building the infrastructure to go further.</p>

  <hr />

  <h3>⚠️ Disclaimer</h3>
  <p>This is a beta release. It's fast. It's wild. It might break. Use at your own risk.</p>

</body>
</html>
