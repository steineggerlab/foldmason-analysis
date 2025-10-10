#!/bin/bash -e

export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max)

# Compute SP/TC/CS scores and extract LDDT/runtime per tool
# ./compute_scores.sh family/ >> scores.tsv

DIR="$1"
FAMILY=$(basename "$DIR")
REF="${DIR}/${FAMILY}_msa.fasta"
AA="${DIR}/sequence.fa"
PDB=$(realpath "${DIR}/pdbs")

# 4 column TSV family, tool, type (sp_fwd/sp_rev/tc/cs), score
compute_score () {
	if [ ! -e "$1" ]; then return; fi
	SCORES="${DIR}/${2}.tcoffee_scores"
	if [ ! -e "$SCORES" ]; then
		SPF=$(t_coffee -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
		SPR=$(t_coffee -other_pg aln_compare -al2 "$REF" -al1 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
		CS=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode column | awk 'NR==3 {print $4}')
		TC=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode tc     | awk 'NR==3 {print $4}')
		printf "%s\t%s\tsp_fwd\t%f\n%s\t%s\tsp_rev\t%f\n%s\t%s\tcs\t%f\n%s\t%s\ttc\t%f\n" \
			"$FAMILY" "$2" "$SPF" \
			"$FAMILY" "$2" "$SPR" \
			"$FAMILY" "$2" "0" \
			"$FAMILY" "$2" "0" > "$SCORES"
	fi
	cat "$SCORES"
}

apply_fn() {
	local func_name=$1
	"$func_name" "${DIR}/foldmason_aa.fa" "foldmason"
	"$func_name" "${DIR}/foldmason_fast_aa.fa" "foldmason_fast"
	"$func_name" "${DIR}/foldmason_refine100_aa.fa" "foldmason_refine100"
	"$func_name" "${DIR}/foldmason_refine100_rosc_aa.fa" "foldmason_refine100_rosc"
	"$func_name" "${DIR}/clustalo.fa" "clustalo"
	"$func_name" "${DIR}/famsa.fa" "famsa"
	"$func_name" "${DIR}/muscle.fa" "muscle"
	"$func_name" "${DIR}/mafft.fa" "mafft"
	"$func_name" "${DIR}/mafft_3di_aa.fa" "mafft_3di"
	"$func_name" "${DIR}/caretta_results/result.fasta" "caretta"
	"$func_name" "${DIR}/mTM_result/result.fasta" "mtmalign"
	"$func_name" "${DIR}/usalign.fa" "usalign"
	"$func_name" "${DIR}/matt.fasta" "matt"
	"$func_name" "${DIR}/mustang.afasta" "mustang"
	"$func_name" "${DIR}/3dcoffee.fa" "3dcoffee"
}


# If the directory has a family_msa.fa, assume it is Homstrad and compute SP/TC/CS
if [ -e "$REF" ]
then
	apply_fn compute_score
	apply_fn compute_nirmsd
fi

# Find all msa2lddt HTML reports and extract LDDT scores
# find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.html" -exec ./extractLDDT.awk {} \;
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name '*.html' -exec ./extractLDDT.awk {} \;

# Get run times for each tool
# find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.time" -exec ./extractTime.awk {} \;
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name '*.time' -exec ./extractTime.awk {} \;
