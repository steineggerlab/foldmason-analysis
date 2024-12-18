#!/usr/bin/awk -f

# Extract amino acid sequence from protein structure PDB file
# ./threetoone.awk file.pdb

BEGIN {
	amino_acid["ALA"] = "A"; amino_acid["CYS"] = "C"; amino_acid["ASP"] = "D"; amino_acid["GLU"] = "E";
	amino_acid["PHE"] = "F"; amino_acid["GLY"] = "G"; amino_acid["HIS"] = "H"; amino_acid["ILE"] = "I";
	amino_acid["LYS"] = "K"; amino_acid["LEU"] = "L"; amino_acid["MET"] = "M"; amino_acid["ASN"] = "N";
	amino_acid["PRO"] = "P"; amino_acid["GLN"] = "Q"; amino_acid["ARG"] = "R"; amino_acid["SER"] = "S";
	amino_acid["THR"] = "T"; amino_acid["VAL"] = "V"; amino_acid["TRP"] = "W"; amino_acid["TYR"] = "Y";

	# Parse file name to use as header, removing paths and .pdb
	name=ARGV[1];
	gsub(".pdb", "", name);
	n = split(name, parts, /\//);
	header = parts[n];
	print ">" header;
}

/^ATOM/ && $3 == "CA" {
	residue = amino_acid[$4];
	if (residue != "") {
		printf "%s", residue;
	}
}

END {
	print ""
}
