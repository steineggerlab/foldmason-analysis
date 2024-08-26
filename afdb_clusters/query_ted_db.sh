#!/bin/bash -e
# ./query_teds.sh ted.sqlite3 afdb.sqlite3 output.tsv

sqlite3 "$1" << EOF
ATTACH DATABASE "$2" AS afdb;
.mode tabs
.output "$3"
.read query_teds.sql
EOF