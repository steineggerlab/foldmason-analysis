#!/usr/bin/awk -f

# Extracts 50x 2-4 domain clusters with 20 members from SQL result

BEGIN {
	OFS="\t"
}

$4 <= 4 {
	count[$4]++;
	if (count[$4] <= 50) {
		split($5, names, ",");

		# Make every cluster representative + 19 others (total 20) 
		final_names[1] = $1;
		num_names = 1;
		for (i = 1; i <= length(names) && num_names < 20; i++) {
			if (names[i] != $1) {
				final_names[++num_names] = names[i];
			}
		}

		# Print row
		printf "%s\t%s\t%s\t%s\t", $1, $2, $3, $4;
		for (i = 1; i <= num_names; i++) {
			printf "%s", final_names[i];
			if (i < num_names) {
				printf ",";
			}
		}
		print "";
	}
}
