#!/bin/bash -e

sqlite3 "$1" << EOF
CREATE TABLE ted_domain(accession text, c text, a text, t text, h text);
.mode tabs
.import "$2" ted_domain
EOF