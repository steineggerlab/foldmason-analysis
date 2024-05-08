#!/bin/bash -e

# families/5_3_exo/foldmason.html

DIR=${1%/*}
FAMILY=${DIR#*/}
FAMILY=${FAMILY##*/}
TOOL=$(basename "$1")
TOOL=${TOOL%\.*}
REF="${DIR}/${FAMILY}_msa.fasta"

SPF=$(t_coffee -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
SPR=$(t_coffee -other_pg aln_compare -al2 "$REF" -al1 "$1" -compare_mode sp     | awk 'NR==3 {print $4}')
CS=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode column | awk 'NR==3 {print $4}')
TC=$(t_coffee  -other_pg aln_compare -al1 "$REF" -al2 "$1" -compare_mode tc     | awk 'NR==3 {print $4}')

# 4 column TSV family, tool, 'lddt', lddt score
printf "%s\t%s\tsp_fwd\t%f\n%s\t%s\tsp_rev\t%f\n%s\t%s\tcs\t%f\n%s\t%s\ttc\t%f\n" \
	"$FAMILY" "$TOOL" "$SPF" \
	"$FAMILY" "$TOOL" "$SPR" \
	"$FAMILY" "$TOOL" "$CS" \
	"$FAMILY" "$TOOL" "$TC"
