#!/usr/bin/awk -f

# Parses time files for tool runs
# Expects FILENAME like families/FAMILY/TOOL.html
# Then reports 4 column TSV family, tool, 'time', lddt score

$0 ~ /^Wall/ {
	n = split(FILENAME, a, "/");
	sub(".time", "", a[n]);
	print a[n-1] "\t" a[n] "\ttime\t" $5;
}
