#!/bin/bash -e

export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max)

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

TOOL_THREADS="${TOOL_THREADS:=1}"
RUN_THREADS="${RUN_THREADS:=1}"
SCORE_THREADS="${SCORE_THREADS:=1}"

# Run all aligners on families in $1
echo "Running aligners"
find $1 -mindepth 1 -maxdepth 1 -type d |\
	THREADS="$TOOL_THREADS" xargs -I{} -P"$RUN_THREADS" ./align_family.sh {}

DIR="$1"
N=$(find $1 -mindepth 1 -maxdepth 1 -type d | wc -l)

check_msas() {
	cnt=$(find "$DIR" -type f -path "*/${1}" | wc -l)
	if [ ! $cnt -eq $N ]; then
		echo "Missing ${1}, ${cnt}/${N}"
		exit 1
	fi
}

echo "Checking all MSAs have been generated"
check_msas "foldmason_aa.fa"
check_msas "foldmason_fast_aa.fa"
check_msas "foldmason_refine100_aa.fa"
check_msas "clustalo.fa"
check_msas "famsa.fa"
check_msas "muscle.fa"
check_msas "mafft.fa"
check_msas "mafft_3di.fa"
check_msas "caretta_results/result.fasta"
check_msas "mTM_result/result.fasta"
check_msas "usalign.fa"
check_msas "matt.fasta"
check_msas "mustang.afasta"
check_msas "3dcoffee.fa"

# echo "Checking all LDDT HTMLs have been generated"
check_msas "foldmason.html"
check_msas "foldmason_fast.html"
check_msas "foldmason_refine100.html"
check_msas "clustalo.html"
check_msas "famsa.html"
check_msas "muscle.html"
check_msas "mafft.html"
check_msas "mafft_3di.html"
check_msas "caretta.html"
check_msas "mtmalign.html"
check_msas "usalign.html"
check_msas "matt.html"
check_msas "mustang.html" #check_msas "3dcoffee.html"

# Get scores per tool
echo "Computing scores"
find $1 -mindepth 1 -maxdepth 1 -type d |\
	xargs -I{} -P"$SCORE_THREADS" ./compute_scores.sh {} |\
       	sort > "$2"
