# AFDB clusters benchmark
## Files
| File | Description |
| --- | --- |
| `final_clusters.tsv` | AlphaFold clusters used in this benchmark (Columns: `rep_accession`, `cath`, `num_members`, `num_domains`, `members`) |
| `filter_ted_domains.awk` | Pre-process TED data for SQL entry (remove AF tags, secondary CATH annotations) |
| `build_ted_db.sh` | Build TED/CATH SQLite3 database from `filter_ted_domains.awk` output |
| `query_ted_db.sh` | Wrapper script for querying TED/CATH database (requires built AFDB SQLite3 database) |
| `query_teds.sql` | Internal SQL script for the cluster query |
| `extract_clusters.awk` | Select clusters from query result |
| `download.sh` | Download PDB files from AlphaFold based on `final_clusters.tsv` |
| `alignments.tar.gz` | Alignments for all tools on AFDB clusters |


## Data preparation
Download and decompress `2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz` and
`5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz` from the [AFDB data dump](https://afdb-cluster.steineggerlab.workers.dev/),
then run the [SQLite3 database build script](https://github.com/steineggerlab/afdb-clusters-web/blob/main/data/build.sh) like:
```sh
wget https://afdb-cluster.steineggerlab.workers.dev/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
wget https://afdb-cluster.steineggerlab.workers.dev/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
pigz -d 2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
pigz -d 5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
./build.sh afdb-clusters.sqlite3 5-allmembers-repId-entryId-cluFlag-taxId.tsv 2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv
```
> Note: this will produce a ~22gb database

Worth running the following commands once this is built:
```sql
PRAGMA optimize;
VACUUM;
```

Download `novel_folds_set.domain_summary.tsv` TED domain data from [Zenodo](https://zenodo.org/records/10848710).
Extract accession/CATH domains from `ted_domain_assignments`.
Only selects the first CATH domain per TED annotation.
```sh
./filter_ted_domains.awk ted_domain_assignments > ted_domains_first_split.tsv
```

Import data into an SQLite database:
```sh
./build_ted_db.sh ted.sqlite3 ted_domains_first_split.tsv
```

Attach the AFDB database and select clusters where members:
* have been clustered by Foldseek (`flag = 2`)
* have at least 2 domains
* match representative CATH domains to T level but are unique at H level
```sh
./query_ted_db.sh ted.sqlite3 afdb.sqlite3 all_clusters.tsv
```

Select clusters for benchmark.
Extracts fifty clusters with 20 structures each of 2, 3 and 4-domain proteins.
```sh
./extract_clusters.awk all_clusters.tsv > final_clusters.tsv
```

Download PDB files
```sh
mkdir afdb_clusters && cd afdb_clusters
./download.sh ../final_clusters.tsv
```

## Run analyses
```sh
../align_families.sh afdb_clusters/ scores.tsv
```