#!/bin/bash -e

# Runs all tools then evaluates LDDT and SoP/TC/CS
#
# ./align_families.sh families/ scores.tsv
#
# Where families is a directory of directories for Homstrad families/AFDB clusters
# Each family directory should be of structure:
# folder/
# 	pdbs/
# 		structureA.pdb
# 		structureB.pdb
# 	folder_msa.fa
# 	folder_aa.fa
#
# Will generate alignment files:
# 	folder/matt/result.fasta
# 	folder/caretta_results/result.fasta
# 	folder/mTM_result/result.fasta
# 	folder/mustang/mustang.afasta
# 	folder/foldmason_aa.fa
# 	folder/clustalo.fa
# 	folder/famsa.fa
# 	folder/mafft.fa
# 	folder/muscle.fa
#
# and HTML reports:
# 	folder/<tool>.html
# 
# scores.tsv will contain LDDT/SoP fwd and rev/TC/CS scores and time in 4 column TSV (family tool scoreType score)
# e.g.
# 	family1	foldmason	lddt 0.6
# 	family1	muscle	sp_fwd	0.3
# 	family1	mafft	tc	0.7

if [ "$#" -ne 2 ]; then
    echo "Error: 2 arguments are required."
    echo "Usage: $0 dataDir/ scores.tsv"
    exit 1
fi

if [ -e "$2" ]; then rm "$2"; fi

THREADS="${THREADS:=1}"

# Run all aligners on families in $1
find $1 -mindepth 1 -maxdepth 1 -type d |\
	THREADS="$THREADS" xargs -I {} ./align_family.sh {} foldmason

# Get scores per tool
for fo in "$1"/*
do
	THREADS="$THREADS" ./compute_scores.sh "$1"
done | sort > "$2"
