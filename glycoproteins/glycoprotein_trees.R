BiocManager::install("YuLab-SMU/treedataverse")

library(ape)
library(Quartet)
library(dplyr)
library(ggtree)
library(ggplot2)
library(svglite)

# Mifsud et al trees
m_e <- read.tree("trees/mifsud_denv1_e_mapped.nw")
m_e1 <- read.tree("trees/mifsud_bvdv1_wsv_e1_mapped.nw")
m_e2 <- read.tree("trees/mifsud_bvdv1_tdav_e2_mapped.nw")
m_e2 <- drop.tip(m_e2, c("HPHV_Hepacivirus_sp."))

# FoldMason trees
f_e <- read.tree("trees/foldmason_denv1_e.nw")  # 3Di tree
# f_e <- read.tree("trees/foldmason_denv1_e_aa.nw")
f_e1 <- read.tree("trees/foldmason_bvdv1_wsv_e1.nw")
f_e2 <- read.tree("trees/foldmason_bvdv1_tdav_e2.nw")

## Remove OKIAV332/WHYC-1 as in their manuscript
f_e <- drop.tip(f_e, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e1 <- drop.tip(f_e1, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e2 <- drop.tip(f_e2, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))

# Split metrics = RF
compute_scores <- function(treeA, treeB) {
  quartets <- QuartetStatus(treeA, treeB)
  splits   <- SplitStatus(treeA, treeB)
  quartet_scores <- SimilarityMetrics(quartets, similarity = TRUE)
  split_scores   <- SimilarityMetrics(splits, similarity = TRUE)
  combined <- rbind(quartet_scores, split_scores)
  row.names(combined) <- c("quartet", "split")
  return(combined) 
}

# Compute all similarity metrics
# Quartet similarity = first row QuartetDivergence
# Robinson-Foulds ≈ second row (split) SymmetricDifference
compute_scores(m_e, f_e)
compute_scores(m_e1, f_e1)
compute_scores(m_e2, f_e2)

# Reproduce tree figures from paper

# a)
# large genome flavi #DEC579
# pesti-like #78777B
# jingmen #82CFF1
# orthoflavi-like #CACACA
# orthoflavi- #217B3A

# b)
# hepaci- #8C1D57
# pegi- #CC6779
# pesti- #44ABA1

clade_colors <- c(
  "largegenome"="#DEC579",
  "pesti-like"="#78777B",
  "jingmen"="#82CFF1",
  "orthoflavi-like"="#CACACA",
  "orthoflavi"="#217B3A",
  "hepaci"="#8C1D57",
  "pegi"="#CC6779",
  "pesti"="#44ABA1"
)

clade_groups <- c(
  "PLPV"="pesti",
  "HPPV"="pegi",
  "HPHV"="hepaci",
  "FJJI"="jingmen",
  "FJMB"="orthoflavi",
  "FJNV"="orthoflavi",
  "FJAF"="orthoflavi",
  "FJFL"="orthoflavi",
  "FJTB"="orthoflavi",
  "FJUN"="orthoflavi",
  "FJIS"="orthoflavi",
  "PLLG"="largegenome",
  "PLUN"="pesti"
)
  
plot_tree <- function(tree, layout="ape") {
  label_clades <- as.data.frame(tree$tip.label, nm="label") %>%
    mutate(clade=clade_groups[sub("_.*$", "", label)])
  p1 <- ggtree(tree, layout=layout, size=0.2) %<+% label_clades +
    geom_tippoint(aes(fill=clade), shape=21, size=1.5, stroke=0.4, color="black") +
    theme_tree() +
    theme(text = element_text(family = "Helvetica")) +
    scale_fill_manual(values = clade_colors) +
    guides(fill="none") +
    geom_treescale(family = "Helvetica", fontsize=2, linesize = 0.2)
  return(p1)
}

plot_tree(m_e) + plot_tree(m_e1) + plot_tree(m_e2) +
  plot_tree(f_e) + plot_tree(f_e1) + plot_tree(f_e2)

ggsave("glycoprotein_trees.pdf", dpi=300)
ggsave("glycoprotein_trees.svg", device=svg, units="mm", width=160, height=100, dpi=300, bg="white")
ggsave("glycoprotein_trees.png", units="mm", width=160, height=100, dpi=300, bg="white")
