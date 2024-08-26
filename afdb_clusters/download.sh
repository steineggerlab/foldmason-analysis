#!/bin/bash -e

# Download PDB files from AlphaFold for each extracted cluster
# Files will be written to <accession>/pdbs/, the expected file structure for align_families.sh

while IFS=$'\t' read -r acc cath nummem numdom members
do
	echo "# Processing $acc"
	mkdir -p "${acc}/pdbs"
	echo "$members" |\
		tr "," "\n" |\
		xargs -I{} echo "https://alphafold.ebi.ac.uk/files/AF-"{}"-F1-model_v4.pdb" |\
		aria2c -s 20 -d "${acc}/pdbs" -i -
done < "$1"
