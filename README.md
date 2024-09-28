<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <h1>Ultidock Project </h1>
</head>
<body>

<p>If you're familiar with docking workflows, this project might come in handy for your research.</p>

<p>All things shorted out now.</p>

<p>To start, download your wget file from your choice of ligand database and put it in your work directory. First run <code>setup.py</code> so that you can create your own workdir as you want. After that, process your ligands using <code>extract.py</code>, then pass the generated files to <code>dock_beta.py</code> -WIP- for docking. This should cover all your docking needs for now. My long-term goal is to automate the entire pipeline and potentially integrate GPU acceleration—though I haven’t found the time for that yet. Be sure to configure the resource allocations in the scripts to align with your system specs for optimal performance. That's my immediate goal, so give it a shot!</p>

<p>Once the docking is complete, my method focuses on finding the ligands with the best affinity, which I define as those in favorable geometric positions and with relatively low binding energy. To achieve this, you need to convert all <code>.pdbqt</code> files into readable data. That’s where <code>output-analyses.py</code> comes in. I opted for CSV format instead of Excel due to the sheer volume of data, but feel free to choose what works best for you. I used Pandas to efficiently sort through all the information.</p>

<p>Since this process demands high-speed random read/write operations, traditional hard drives won’t cut it. You’ll need a fast NVMe SSD, or even better, an Intel Optane drive. <strong>INTEL, ARE YOU LISTENING?</strong></p>

<h3><strong>DISCLAIMER:</strong> This is an experimental project. Use at your own risk.</h3>

<p>For context, I ran simulations with all the ligands from the <code>wget</code> file, which took 3 days. After processing, I ended up with over 1.2 million ligand files, taking up around 80GB of storage. My setup is modest—a Ryzen 5 3600X with 24GB of RAM. If you’ve got access to a serious server with many cores and fast NVMe storage (Optane, perhaps), please reach out. I’d love to run simulations for a wider range of molecules.</p>

<p>For now, I’ve uploaded my attempt at docking with 4H10. I identified a few promising candidates, but keep in mind—I’m neither a molecular physicist nor a bioinformatician. Just a physicist navigating toward more advanced simulations.</p>

</body>
</html>
