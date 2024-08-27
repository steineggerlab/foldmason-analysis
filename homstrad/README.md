# Homstrad benchmark
## Files
| File | Description |
| --- | --- |
| `clean_homstrad.py` | Clean Homstrad family folders, generate individual PDB files, AA FASTA files, Homstrad MSA |
| `download_homstrad.sh` | Downloads latest Homstrad database from FTP, unpacks and fixes igV family |
| `alignments.tar.gz` | Alignments from all tools |

## Data preparation
Run `download_homstrad.sh` to download latest Homstrad release and prepare it for analysis.
Generates directories `homstrad_db` (raw) and `homstrad_clean` (processed).

## Run analyses
Run `./align_families.sh homstrad_clean homstrad_scores.tsv` to run the full suite of
tools on the Homstrad database, and save all scores and runtimes to `homstrad_scores.tsv`.
