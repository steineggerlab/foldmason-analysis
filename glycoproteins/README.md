# foldmason-analysis

| File | Description |
| --- | --- |
| extract_sequences.sh | Extract sequences (+-500aa) from Mifsud et. al. for re-prediction |
| glycoprotein_trees.R | Compute tree similarities and generate visualisations |
| trees/ | Mifsud et. al. and FoldMason tree files |

# Instructions
## Prerequisites
Ensure you have downloaded the Mifsud et. al. dataset from [Zenodo](https://zenodo.org/records/10616318)
and it is unpacked in your working directory. It should resemble:

```
folder/
    Mifsud_2024_data/
        structure guided glycoprotein phylogenies/
        reference structures/
        tables/
        trees/
        structures/
        sequences/
        alignments/
    extract_sequences.sh
```

## Process
1. Run `extract_sequences.sh` to extract sequences for ColabFold prediction. This will
   extract +-500aa of the first/last aligned blocks in the Mifsud et. al. MSAs.
