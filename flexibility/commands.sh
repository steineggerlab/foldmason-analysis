# Run inside folder like
# family/
#   all/
#     1.pdb
#     2.pdb
#     3.pdb
#   pair/
#     1.pdb
#     2.pdb

# FoldMason
# Will reuse db for msa2lddt
foldmason easy-msa all fm_msa tmp
foldmason easy-msa pair fm_pair tmp_pair
foldmason msa2lddt tmp/latest/structures fm_msa_aa.fa
foldmason msa2lddt tmp/latest/structures fm_pair_aa.fa

# Matt
Matt $(find pair -type f -name '*.pdb') -o matt_pair -s0 -b1
Matt $(find all -type f -name '*.pdb') -o matt_msa -s0 -b1
sed -i 's/:.*$//g' matt_msa_bent.fasta
sed -i 's/:.*$//g' matt_pair_bent.fasta
foldmason msa2lddt tmp/latest/structures matt_pair_bent.fasta
foldmason msa2lddt tmp/latest/structures matt_msa_bent.fasta

# Caretta
caretta-cli pair -o caretta_pair
caretta-cli all -o caretta_msa
sed -i 's/\.pdb//g' caretta_pair/result.fasta
sed -i 's/\.pdb//g' caretta_msa/result.fasta
foldmason msa2lddt tmp/latest/structures caretta_pair/result.fasta
foldmason msa2lddt tmp/latest/structures caretta_msa/result.fasta

# mTM-align
mTM-align -i <(ls all) -o mtm_msa
mTM-align -i <(ls pair) -o mtm_pair
sed -i 's/\.pdb//g' mtm_msa/result.fasta
sed -i 's/\.pdb//g' mtm_pair/result.fasta
foldmason msa2lddt tmp/latest/structures mtm_msa/result.fasta
foldmason msa2lddt tmp/latest/structures mtm_pair/result.fasta

# MUSTANG
mustang -p all -o mustang_msa -F fasta
mustang -p pair -o mustang_pair -F fasta
sed -i 's/\.pdb//g' mustang_msa.afasta
sed -i 's/\.pdb//g' mustang_pair.afasta
foldmason msa2lddt tmp/latest/structures mustang_msa.afasta
foldmason msa2lddt tmp/latest/structures mustang_pair.afasta

# US-align
USalign pair/*.pdb -outfmt 1 > usalign_pair.fasta
USalign all <(find all -type f -name '*.pdb') -outfmt 1 > usalign_msa.fasta
sed -i -e '/^[#$]/d; /^$/d; s/\.pdb.*$//g; s/>.*\//>/' usalign_pair.fasta
sed -i -e '/^[#$]/d; /^$/d; s/\.pdb.*$//g; s/>.*\//>/' usalign_msa.fasta
foldmason msa2lddt tmp/latest/structures usalign_msa.afasta
foldmason msa2lddt tmp/latest/structures usalign_pair.afasta