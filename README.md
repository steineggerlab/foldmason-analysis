# foldmason-analysis

| File | Description |
| --- | --- |
| clean_homstrad.py | Clean Homstrad family folders, generate individual PDB files, AA FASTA files, Homstrad MSA |
| download_homstrad.sh | Downloads latest Homstrad database from FTP, unpacks and fixes igV family |
| compute_scores.sh | Compute SoP (forward and backward), TC and CS scores (Homstrad), extract LDDT and runtime |
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

Also that the directory structure of the dataset resembles:

```
folder/
    family1/
        pdbs/
            structure1.pdb
            structure2.pdb
            structure3.pdb
            ...
        family1_msa.fa (Homstrad reference alignment)
        family1_aa.fa
    family2/
        ...
    ...
```

## Process
1. Run `download_homstrad.sh` to download latest Homstrad release and prepare it for analysis.
   Generates directories `homstrad_db` (raw) and `homstrad_clean` (processed).
2. Run `./align_families.sh homstrad_clean homstrad_scores.tsv` to run the full suite of
   tools on the Homstrad database, and save all scores and runtimes to `homstrad_scores.tsv`.
