#!/bin/bash -ex
# Compute all-vs-all alignments of protein families using MMseqs2 and Foldseek.
# Used to generate data for Figure S5 (AFDB cluster identity distributions)

fake_pref() {
	QDB="$1"
	TDB="$2"
	RES="$3"
	ln -s "${TDB}.index" "${RES}"
	INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
	awk -v size=$INDEX_SIZE '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
	awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

family=$(basename "$1")
db=$(readlink -f "${1}/foldmason/latest/structures")
prefdb="${1}/fakepref"
mm_alndb="${1}/mm_allvsalldb"
mm_res="${1}/mm_allvsall.m8"
fs_alndb="${1}/fs_allvsalldb"
fs_res="${1}/fs_allvsall.m8"
output="${1}/allvsall.csv"

if [ -L "$prefdb" ]; then
	rm "$prefdb"
fi
fake_pref "$db" "$db" "$prefdb"
mmseqs align "$db" "$db" "$prefdb" "$mm_alndb" --threads 1 -e inf > /dev/null 2>&1
mmseqs convertalis "$db" "$db" "$mm_alndb" "$mm_res" \
	--threads 1 \
	--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov > /dev/null 2>&1
foldseek structurealign "$db" "$db" "$prefdb" "$fs_alndb" --threads 1 -e inf > /dev/null 2>&1
foldseek convertalis "$db" "$db" "$fs_alndb" "$fs_res" \
	--threads 1 \
	--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov #> /dev/null 2>&1

awk -v fam="$family" -F'\t' '$3 != 1.0 { print fam",mmseqs,"$1","$2","$3","$13","$14 }'   "$mm_res" > "$output"
awk -v fam="$family" -F'\t' '$3 != 1.0 { print fam",foldseek,"$1","$2","$3","$13","$14 }' "$fs_res" >> "$output"

# Run on all families, then:
# find families_new/ -mindepth 1 -maxdepth 1 -type d | xargs -I{} -P32 ./identity.sh {}
# find families_new/ -mindepth 2 -maxdepth 2 -name 'allvsall.csv' -exec cat {} + > identities.csv
