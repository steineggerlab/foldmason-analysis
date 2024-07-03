#!/bin/bash -e

# From Foldseek Jupyter notebook
# Download latest Homstrad, unpack and fix chain assignments in igV family
# Decompressed homstrad release in homstrad_db/, cleaned in homstrad_clean/

if [ ! -d homstrad_db ]; then
	mkdir homstrad_db
	curl https://homstrad.mizuguchilab.org/homstrad/data/homstrad_with_PDB_2024_May_1.tar.gz | tar -xz -C homstrad_db
	
	# Fix the igV PDB, because the chain order is wrong (corrects only the first two)
	mv homstrad_db/igV/igV-sup.pdb homstrad_db/igV/igV-sup.pdb.orig
	
	echo \
'REMARK The domains in this file are:
REMARK    1tvdb   chain    H
REMARK    1qfpa   chain    G
REMARK    1b88a   chain    A
REMARK    1cd8    chain    B
REMARK    3cd4    chain    I
REMARK    1neu    chain    F
REMARK    1hnga   chain    E
REMARK    1cid1   chain    C
REMARK    1hnf    chain    D' > homstrad_db/igV/igV-sup.pdb
	cat homstrad_db/igV/igV-sup.pdb.orig >> homstrad_db/igV/igV-sup.pdb
fi

if [ ! -d homstrad_clean ]; then
	find homstrad_db -mindepth 1 -maxdepth 1 -type d | xargs -I{} -P4 bash -c 'python3 clean_homstrad.py "$1" "${1/homstrad_db/homstrad_clean}"' - '{}'
fi
