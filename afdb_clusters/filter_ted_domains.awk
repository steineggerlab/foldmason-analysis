#!/usr/bin/awk -f

/^AF/ {
    acc=$1;
    cath=$6;
    sub(/AF-/, "", acc);
    sub(/-F1-model_v4.*/, "", acc);
    sub(/,.*$/, "", cath);
    gsub(/\./, "\t", cath);
    print acc,cath;
}