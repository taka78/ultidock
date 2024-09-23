If you know what you're doing, this could be handy for your research.

You're gonna need the Autodock Vina sources on your computer, and yeah, you'll have to deal with the directory mess I created. But don't worry, it's not that bad.

First, prep your ligands with extract.py, then run dock_beta.py with the output files. That should handle all your docking needs. My ultimate goal is to automate the whole process and eventually add GPU acceleration, but I haven't found the time for that yet. Make sure to configure the resources in all the scripts to match your hardware setup for maximum performance. That's basically my first target, so go ahead and DO IT!

THIS IS NOT A STABLE PROJECT. USE IT AT YOUR OWN RISK.

BTW, for comparison, I ran simulations with all the ligands from the wget file over 3 days. After separating them, there were over 1.2 million ligand files, totaling around 80GB in size. My computer is just a basic machine with Ryzen 5 3600X with 24GB of RAM. If you have a real server with lots of cores and NVMe SSDs —maybe even Optane— please let me know so I can run simulations for all known molecules.

I've uploaded my 4H10 docking attempt. I found some possible suspects, but I'm neither a molecular physicist nor a bioinformatician (?). I'm just a simple physicist, trying to find my way to ultimate simulations.
