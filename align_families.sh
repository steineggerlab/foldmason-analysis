#!/bin/bash -e

# Runs all tools then evaluates LDDT and SoP/TC/CS
#
# ./run_homstrad.sh families/ scores.tsv times.tsv
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
# scores.tsv will contain LDDT/SoP/TC/CS scores in 4 column TSV (family tool scoreType score)
# e.g.
# 	family1 foldmason lddt 0.6
# 	family1 muscle sp_fwd 0.3
# 	family1 mafft tc 0.7

if [ "$#" -ne 3 ]; then
    echo "Error: 3 arguments are required."
    echo "Usage: $0 dataDir/ scores.tsv times.tsv"
    exit 1
fi

if [ -e "$2" ]; then rm "$2"; fi
if [ -e "$3" ]; then rm "$3"; fi

THREADS="${THREADS:=1}"

# Run all aligners on families in $1
find $1 -mindepth 1 -maxdepth 1 -type d |\
	xargs -I {} THREADS="$THREADS" ./run_family.sh {}

# Score everything; take LDDT score from HTML reports, compute SP/TC/CS with T-coffee
for fo in "$1"/*
do
	if [ ! -d "$fo" ]; then
		continue
	fi
	find "$fo" -maxdepth 1 -name '*.html' | xargs -I{} -P "$THREADS" extractLDDT.awk {} >> "$2"

	# If the directory has a family_msa.fa, assume it is Homstrad and compute SP/TC/CS
	family=$(basename "$fo")
	if [ -e "${fo}/${basename}_msa.fa" ]; then
		./compute_spcstc.sh "${fo}/foldmason_aa.fa" >> "$2"
		./compute_spcstc.sh "${fo}/clustalo.fa" >> "$2"
		./compute_spcstc.sh "${fo}/famsa.fa" >> "$2"
		./compute_spcstc.sh "${fo}/muscle.fa" >> "$2"
		./compute_spcstc.sh "${fo}/mafft.fa" >> "$2"
		./compute_spcstc.sh "${fo}/caretta_results/result.fasta" >> "$2"
		./compute_spcstc.sh "${fo}/mTM_result/result.fasta" >> "$2"
		./compute_spcstc.sh "${fo}/matt/matt.fasta" >> "$2"
		./compute_spcstc.sh "${fo}/mustang/mustang.afasta" >> "$2"
	fi

	# Get run times for each tool
	find "$fo" -maxdepth 1 -type f -name "*.time" |\
		xargs -I{} -P"$THREADS" awk '$0 ~ /^Wall/ {split(FILENAME, a, "/"); sub(".time", "", a[3]); print a[2]"\t"a[3]"\t"$5}' {} >> $3
done
