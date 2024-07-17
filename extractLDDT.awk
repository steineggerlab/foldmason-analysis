#!/usr/bin/awk -f

# Parses msa2lddt result HTML files for LDDT scores
# Expects FILENAME like families/FAMILY/TOOL.html
# Then reports 4 column TSV family, tool, 'lddt', lddt score

# GNU awk
# {
# 	if (match($0, /msaLDDT":([0-9]*\.[0-9]*)[},]/, m)) {
# 		n = split(FILENAME, a, "/");
# 		sub(/\.html/, "", a[n]);
# 		print a[n-1] "\t" a[n] "\tlddt\t" m[1];
# 	}
# }

# POSIX
BEGIN {
	infile=ARGV[1];
	match(infile, /\/[A-Z0-9a-z_]+\.html/);
	tool=substr(infile, RSTART+1, RLENGTH-6);  # remove / and .html
	family=infile
	sub(/\/[^\/]+$/, "", family);  # remove filename
	sub(/.*\//, "", family);       # remove path until last dir
}

{
	# LDDT score
	match($0, /msaLDDT":[0-9]+\.[0-9]+/);
	if (RLENGTH > 0) {
		printf "%s\t%s\tlddt\t%s\n", family,tool,substr($0, RSTART+9, RLENGTH-9)
	}

	# Total columns contributing to score
	match($0, /\"scores\"\: \[.+\]/);
	if (RLENGTH > 0) {
		count = 0;
		array = substr($0, RSTART+11, RLENGTH-12);  # remove "scores": [ and ]
		n = split(array, elements, ",");
		for (i = 1; i <= n; i++) {
			if (elements[i] != "-1") {
				count++;
			}
		}
		printf "%s\t%s\tlddt_columns\t%s/%s\n", family,tool,count,n
	}
}
