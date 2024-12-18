#!/bin/bash -e

# Download PDB files from AlphaFold for each extracted cluster
# Files will be written to <accession>/pdbs/, the expected file structure for align_families.sh
# Will also produce <accession>/sequence.fa by extracting AA sequences from downlaoded PDB files

dir=$(dirname $(readlink -f "$0"))

while IFS=$'\t' read -r acc cath nummem numdom members
do
	echo "# Processing $acc"
	pdbs="${acc}/pdbs"
	mkdir -p "$pdbs"
	echo "$members" |\
		tr "," "\n" |\
		xargs -I{} echo "https://alphafold.ebi.ac.uk/files/AF-"{}"-F1-model_v4.pdb" |\
		aria2c -s 20 -d "${acc}/pdbs" -i -
	find "$pdbs" -type f -name "*.pdb" -exec "${dir}/extract_sequences.awk" {} \; > "${acc}/sequence.fa"
done < "$1"
