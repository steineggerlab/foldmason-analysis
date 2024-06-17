#!/bin/bash -e
# set -x

# ./extract_sequences.sh
# Extracts sequences from Mifsud et al for ColabFold re-prediction
# Takes +-500aa around the aligned block ranges from MSAs

# Proteome sequence file
PROTEINS="../Mifsud_2024_Data/sequences/full_genome_sequences/flaviviridae_protein_seqs.fasta"

# Manually curated mapping file to map sequence headers <=> proteome sequences
# e.g. PLLG_k141_32700_flag1_multi37.0000_len21298_SRR9960026_Phenacoccus_solenopsis_CottonMealybug
#      => PLLG_Phenacoccus_solenopsis_associated_flavi-like_virus
MAPPING="../mifsud_id_mapping.csv"

# Input files
# 251, 191 and 168 sequences respectively
declare -a FILES=(
	"../Mifsud_2024_Data/structure guided glycoprotein phylogenies/E glycoprotein/DENV1_E_reference_foldseek_aligned_20231128.fas"
	"../Mifsud_2024_Data/structure guided glycoprotein phylogenies/E1 glycoprotein/BVDV1_WSV_E1_reference_foldseek_aligned_20231128.fas"
	"../Mifsud_2024_Data/structure guided glycoprotein phylogenies/E2 glycoprotein/BVDV1_TDAV_E2_reference_foldseek_aligned_20231128.fas"
)

# Get sequence lengths
LENGTHS="glyco_lengths.tsv"
awk '/^>/ {if (seq) print header, length(seq); header = substr($0, 2); seq=""} !/^>/ {seq = seq $0} END {print header"\t"length(seq)}' "$PROTEINS" |\
	sort > "$LENGTHS"

# For each MSA:
# 	1. Map non-matching sequence headers to those in proteome file
# 	2. Extract aligned blocks from FASTA header, join to sequence length from proteome
# 	3. Take +-500aa around the block range (min 0, max sequence length)
# 	4. Extract sequence from proteome file
for msa in "${FILES[@]}"; do
	name=$(basename "$msa")
	code="${name%%_reference*}"
	awk -F',' '{ print "s/"$1"/"$2"/g" }' "$MAPPING" | sed -f - "$msa" |\
		awk '/^>/ { sub(/>/, "", $0); sub(/\|.*?\|/, "\t", $0); sub(/\r/, "", $0); print; }' | sort |\
		join -1 1 - "$LENGTHS" > "${code}_aligned_blocks.tsv"
	awk '{ split($2, blocks, "_"); start=blocks[1]*100-500; end=(blocks[length(blocks)])*100+300+500; print $1, (start < 0) ? 0 : start, (end > $3) ? $3 : end }' "${code}_aligned_blocks.tsv" > "${code}_ranges.tsv"
	awk 'FNR==NR && /^>/ { header=substr($0,2); seq[header]=""; next; } FNR==NR { seq[header]=seq[header] $0; next; } /^#/ {next;} {if ($1 in seq) {subseq=substr(seq[$1], $2, $3-$2+1); print ">"$1"_"$2"_"$3"\n"subseq} }' "$PROTEINS" "${code}_ranges.tsv" > "${code}_sequences.fa"
done
