# foldmason-analysis

| File | Description |
| --- | --- |
| clean_homstrad.py | Clean Homstrad family folders, generate individual PDB files, AA FASTA files, Homstrad MSA |
| compute_spcstc.sh | Compute SoP (forward and backward), TC and CS scores for a Homstrad family |
| download_homstrad.sh | Downloads latest Homstrad database from FTP, unpacks and fixes igV family |
| extractLDDT.awk | Extract LDDT scores from msa2lddt HTML report |
| align_family.sh | Run tools on a given folder and generate msa2lddt reports on each resulting MSA |
| align_families.sh | Run tools on all subdirectories of given folder and compute LDDT scores (+SP/TC/CS if Homstrad MSA found) |

# Instructions
## Prerequisites
Ensure the following tools are available on system `$PATH`:

```
foldmason
caretta-cli
muscle (version 5)
famsa (version 2)
mafft (linsi mode)
clustalo
Matt
mustang
mTM-align
t_coffee (for computing SoP/TC/CS scores)
```

## Process
1. Run `download_homstrad.sh` to download latest Homstrad release and prepare it for analysis.
   Generates directories `homstrad_db` (raw) and `homstrad_clean` (processed).
2. Run `./run_homstrad.sh homstrad_clean homstrad_scores.tsv` to run the full suite of
   tools on the Homstrad database and save all scores to `homstrad_scores.tsv`.
