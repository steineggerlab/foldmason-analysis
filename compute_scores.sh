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
	SPF=$(t_coffee -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
	SPR=$(t_coffee -other_pg aln_compare -al2 "$REF" -al1 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
	CS=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode column | awk 'NR==3 {print $4}')
	TC=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode tc     | awk 'NR==3 {print $4}')
	printf "%s\t%s\tsp_fwd\t%f\n%s\t%s\tsp_rev\t%f\n%s\t%s\tcs\t%f\n%s\t%s\ttc\t%f\n" \
		"$FAMILY" "$2" "$SPF" \
		"$FAMILY" "$2" "$SPR" \
		"$FAMILY" "$2" "$CS" \
		"$FAMILY" "$2" "$TC"
}

compute_nirmsd() {
	if [ ! -e "$1" ]; then return; fi
	TEMPLATE="${DIR}/tcoffee.template"
	LOG="${DIR}/${2}_nirmsd.log"
	if [ ! -e "$TEMPLATE" ]; then
		awk -v fo=$(realpath "$PDB") \
			'/^>/ { header=$1; sub(/^>/, "", header); print $1" _P_ "fo"/"header".pdb" }' \
			"$AA" > "$TEMPLATE"
	fi
	if [ ! -e "$LOG" ]; then
		t_coffee -other_pg irmsd "$1" -template_file "$TEMPLATE" -quiet /dev/null > "$LOG"
	fi
	if grep -q "ERROR:" "$LOG" || grep -q "FATAL:T-COFFEE:" "$LOG"; then
		t_coffee -other_pg irmsd "$1" -template_file "$TEMPLATE" -quiet /dev/null > "$LOG"
	fi
	awk -v fam="$FAMILY" -v tool="$2" '
		/TOTAL\s*APDB:/   {print fam "\t" tool "\tapdb\t" $3}
		/TOTAL\s*iRMSD:/  {print fam "\t" tool "\tirmsd\t" $3}
		/TOTAL\s*NiRMSD:/ {print fam "\t" tool "\tnirmsd\t" $3}
	' "$LOG"
}

apply_fn() {
	local func_name=$1
	"$func_name" "${DIR}/foldmason_aa.fa" "foldmason"
	"$func_name" "${DIR}/foldmason_refine100_aa.fa" "foldmason_refine100"
	"$func_name" "${DIR}/clustalo.fa" "clustalo"
	"$func_name" "${DIR}/famsa.fa" "famsa"
	"$func_name" "${DIR}/muscle.fa" "muscle"
	"$func_name" "${DIR}/mafft.fa" "mafft"
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
       compute_nirmsd "$REF" "homstrad" 
fi

# Compute NiRMSD with T-Coffee
apply_fn compute_nirmsd

# Find all msa2lddt HTML reports and extract LDDT scores
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.html" -exec ./extractLDDT.awk {} \;

# Get run times for each tool
find "$DIR" -mindepth 1 -maxdepth 1 -type f -name "*.time" -exec ./extractTime.awk {} \;
