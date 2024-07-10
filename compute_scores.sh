#!/bin/bash -e

# Compute SP/TC/CS scores and extract LDDT/runtime per tool
# ./compute_scores.sh family/ >> scores.tsv

DIR="$1"
FAMILY=$(basename "$DIR")
REF="${DIR}/${FAMILY}_msa.fasta"

# 4 column TSV family, tool, type (sp_fwd/sp_rev/tc/cs), score
compute_score () {
	if [ ! -e "$1" ]; then return; fi
	TOOL=$(basename "$1" | sed 's/_aa\.fa//g')
	TOOL=${TOOL%\.*}
	SPF=$(t_coffee -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
	SPR=$(t_coffee -other_pg aln_compare -al2 "$REF" -al1 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
	CS=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode column | awk 'NR==3 {print $4}')
	TC=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode tc     | awk 'NR==3 {print $4}')
	printf "%s\t%s\tsp_fwd\t%f\n%s\t%s\tsp_rev\t%f\n%s\t%s\tcs\t%f\n%s\t%s\ttc\t%f\n" \
		"$FAMILY" "$TOOL" "$SPF" \
		"$FAMILY" "$TOOL" "$SPR" \
		"$FAMILY" "$TOOL" "$CS" \
		"$FAMILY" "$TOOL" "$TC"
}

# If the directory has a family_msa.fa, assume it is Homstrad and compute SP/TC/CS
if [ -e "$REF" ]
then
	compute_score "${DIR}/foldmason_aa.fa"
	compute_score "${DIR}/foldmason_refine1000_aa.fa"
	compute_score "${DIR}/clustalo.fa"
	compute_score "${DIR}/famsa.fa"
	compute_score "${DIR}/muscle.fa"
	compute_score "${DIR}/mafft.fa"
	compute_score "${DIR}/caretta_results/result.fasta"
	compute_score "${DIR}/mTM_result/result.fasta"
	compute_score "${DIR}/matt/matt.fasta"
	compute_score "${DIR}/mustang/mustang.afasta"
fi

# Find all msa2lddt HTML reports and extract LDDT scores
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.html" -exec ./extractLDDT.awk {} \;

# Get core (ungapped) LDDT scores
# find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.html" -exec ./getCoreLDDT.py {} \;

# Get run times for each tool
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.time" -exec ./extractTime.awk {} \;
