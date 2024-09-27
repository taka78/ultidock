#!/bin/bash

# Get a list of all the files to process
FILES="/mnt/g/ligands-prep/*"


# Start a parallel job for each file
for f in $FILES
do
    echo "Processing $f"
    ./vina_split --input /mnt/g/ligands-prep/$f.pdbqt
    rm /mnt/g/ligands-prep/$f.pdbqt

