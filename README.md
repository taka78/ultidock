<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
</head>
<body>

  <h1>Ultidock – GROMACS Integration Branch (<code>gmx-dev</code>)</h1>

  <p>
    <strong>Ultidock</strong> is a high-throughput docking and simulation workflow. This branch focuses on
    reliable automation and future GROMACS integration, with particular attention to modern NVIDIA GPUs.
  </p>

  <div class="note">
    <strong>Key points:</strong>
    <ul>
      <li><strong>CUDA Toolkit ≥ 12.8</strong> is required for modern GPUs (e.g., RTX 40/50 series, datacenter parts).</li>
      <li><strong>All setup</strong> (GPU/driver checks, architecture targets, compilation, directory creation) is handled by <code>run.py</code>.</li>
      <li><strong>User actions:</strong> place your receptor file in <code>MACRO_MOL_DIR/</code> and replace <code>ligands.wget</code> with your ZINC22 (or preferred database) list.</li>
    </ul>
  </div>

  <hr/>

  <h2>1. Features</h2>
  <ul>
    <li><strong>Automatic build &amp; environment setup:</strong> <code>run.py</code> detects NVIDIA GPUs and chooses appropriate compilation targets. It compiles AutoDock-GPU and prepares required directories.</li>
    <li><strong>Deterministic output management:</strong> docking outputs are written directly into <code>DOCKING_DIR/</code> with unique, parse-ready names.</li>
    <li><strong>Thread-safety by design:</strong> GPU runs are serialized per device to avoid CUDA context races; CPU preparation remains parallel.</li>
    <li><strong>Extensibility:</strong> planned integration with GROMACS and optional visualization layers (OpenGL/Web) without impacting the core pipeline.</li>
  </ul>

  <h2>2. Required Software &amp; Libraries</h2>
  <h3>2.1 Core</h3>
  <ul>
    <li><strong>NVIDIA Driver</strong> (sufficiently recent for CUDA ≥ 12.8)</li>
    <li><strong>CUDA Toolkit 12.8</strong> (or newer)</li>
    <li><strong>GNU Toolchain:</strong> <code>gcc</code>/<code>g++</code>, <code>make</code></li>
    <li><strong>Python 3.8+</strong> &nbsp;packages: <code>numpy</code>, <code>biopython</code>, <code>pandas</code>, <code>tqdm</code></li>
    <li><strong>AutoDock-GPU</strong> (compiled by <code>run.py</code> into <code>AUTODOCK_GPU_DIR/bin/</code>)</li>
    <li><strong>AutoDock Vina</strong> (CPU fallback, optional)</li>
    <li><strong>AutoGrid</strong> (compiled by the accompanying script if needed)</li>
  </ul>

  <h3>2.2 Optional / Future</h3>
  <ul>
    <li><strong>GROMACS</strong> (GPU build recommended) for post-docking MD and free-energy analysis</li>
    <li><strong>OpenGL</strong> development headers for native visualization</li>
    <li><strong>MPI</strong> if deploying to clusters in later stages</li>
  </ul>

  <hr/>

  <h2>3. Directory Layout</h2>
  <pre><code>workdir/
 ├─ data-analyses
 ├─ docking
 └├─ run.py
  ├─ ligands.wget                # replace with your list (e.g., ZINC22)
  ├─ MACRO_MOL_DIR/              # place your receptor .pdbqt here
  ├─ LIGANDS_DIR/                # auto-created; ligands fetched/organized
  ├─ DOCKING_DIR/                # docking outputs (_out.pdbqt) land here
  ├─ AUTODOCK_GPU_DIR/           # AutoDock-GPU sources and built binaries
  └─ RESULTS_DIR/                # summary data, logs, analysis artifacts
</code></pre>

  <hr/>

  <h2>4. Quick Start</h2>
  <ol>
    <li>Ensure the NVIDIA driver and CUDA Toolkit <strong>12.8+</strong> are installed.</li>
    <li>Copy your receptor file (e.g., <code>protein.pdbqt</code>) into <code>MACRO_MOL_DIR/</code>.</li>
    <li>Replace <code>ligands.wget</code> with a list from <strong>ZINC22</strong> (or your preferred source).</li>
    <li>Run:
      <pre><code class="language-bash">/usr/bin/python3 run.py</code></pre>
    </li>
  </ol>

  <p><strong>What happens:</strong> <code>run.py</code> will validate the environment, compile AutoDock-GPU (selecting target architectures automatically), prepare directories, download/process ligands, and run docking.</p>

  <hr/>

  <h2>5. Usage Notes</h2>
  <ul>
    <li><strong>GPU-CPU Fallback:</strong> If you don't have a supported GPU, parallelized CPU docking will start with Vina.</li>
    <li><strong>GPU scheduling:</strong> docking tasks are dispatched one-at-a-time per GPU to ensure stability. CPU-side preparation may remain parallel.</li>
    <li><strong>Output format:</strong> AutoDock-GPU writes <code>&lt;basename&gt;_out.pdbqt</code> under <code>DOCKING_DIR/</code>. Filenames include a unique suffix to avoid collisions in batch runs.</li>
    <li><strong>Repeatability:</strong> compilation and directory creation are idempotent; rerunning <code>run.py</code> is safe.</li>
  </ul>

  <hr/>

  <h2>6. Troubleshooting</h2>
  <ul>
    <li><strong>Toolkit mismatch:</strong> if compilation fails with “sm_&lt;arch&gt; not supported”, upgrade to CUDA ≥ 12.8 and re-run <code>run.py</code>.</li>
    <li><strong>Parallelism errors:</strong> if you observe CUDA initialization assertions, ensure only one docking runs per GPU at any moment. <code>run.py</code> enforces this by design.</li>
    <li><strong>Paths:</strong> use absolute paths if invoking binaries directly. The pipeline already passes explicit file paths to avoid CWD-related issues.</li>
  </ul>

  <hr/>

  <h2>7. Future Updates</h2>
  <ul>
    <li><strong>Multi-GPU scaling:</strong> per-GPU workers with automatic enumeration; multi-node options under evaluation.</li>
    <li><strong>GROMACS integration:</strong> automated post-docking MD setup and free-energy estimation.</li>
    <li><strong>Visualization:</strong> optional OpenGL/Desktop viewer and WebGL dashboard for pose inspection.</li>
    <li><strong>Mac/CPU paths:</strong> streamlined CPU fallback for environments without NVIDIA GPUs.</li>
  </ul>

  <hr/>

  <h2>8. Change Summary (this branch)</h2>
  <ul>
    <li>Auto-detection of GPU architecture and compatible <code>TARGETS</code> during build.</li>
    <li>Stable output routing to <code>DOCKING_DIR/</code> and unique file naming.</li>
    <li>Thread-safety and resource control for reliable GPU execution.</li>
    <li>Documentation aligned to CUDA 12.8+ and automated setup via <code>run.py</code>.</li>
  </ul>
  <h2>9. Citation</h2>
  
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

<h2>10. Acknowledgements</h2>
  <p>
  Ultidock relies on the robust and widely used <a href="http://vina.scripps.edu">AutoDock Vina</a> software for molecular docking. If you use Ultidock, please also cite the original Vina publication:
  </p>

  <blockquote>
  Trott, O., & Olson, A. J. (2010). <em>AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.</em> Journal of Computational Chemistry, 31(2), 455–461.  
  <a href="https://doi.org/10.1002/jcc.21334">https://doi.org/10.1002/jcc.21334</a>
  </blockquote>

  <blockquote>
  <p>
  Diogo Santos-Martins, Leonardo Solis-Vasquez, Andreas F Tillack, Michel F Sanner, Andreas Koch, and Stefano Forli (2021)
  Accelerating AutoDock4 with GPUs and Gradient-Based Local Search
  Journal of Chemical Theory and Computation 2021 17 (2), 1060-1073  <p>
  https://doi.org/10.1021/acs.jctc.0c01006
  <hr />

<p align=\"center\">If you find Ultidock useful, please star the repository!</p>

</body>
</html>

</body>
</html>