#!/usr/bin/awk -f

# Parses msa2lddt result HTML files for LDDT scores
# Expects FILENAME like families/FAMILY/TOOL.html
# Then reports 4 column TSV family, tool, 'lddt', lddt score

{
	if (match($0, /msaLDDT":([0-9]*\.[0-9]*)[},]/, m)) {
		split(FILENAME, a, "/");
		sub(/\.html/, "", a[3]);
		print a[2] "\t" a[3] "\tlddt\t" m[1]
	}
}
