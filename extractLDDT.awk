#!/usr/bin/awk -f

# Parses msa2lddt result HTML files for LDDT scores
# Expects FILENAME like families/FAMILY/TOOL.html
# Then reports 4 column TSV family, tool, 'lddt', lddt score

{
	if (match($0, /msaLDDT":([0-9]*\.[0-9]*)[},]/, m)) {
		n = split(FILENAME, a, "/");
		sub(/\.html/, "", a[n]);
		print a[n-1] "\t" a[n] "\tlddt\t" m[1]
	}
}
