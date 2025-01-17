# For TreeIO/ggtree
# BiocManager::install("YuLab-SMU/treedataverse")

library(Quartet)
library(ape)
library(dplyr)
library(tidyr)
library(ggtree)
library(ggplot2)
library(svglite)
library(patchwork)
library(grid)

BASEDIR="./glycoproteins/"

# Flavi-Jingmen Clade
# FJTB = Flavi-Jingmen Tick-Borne
# FJMB = Flavi-Jingmen Mosquito-Borne
# FJNV = Flavi-Jingmen No Known Vector
# FJIS = Flavi-Jingmen Insect Only
# FJAF = Flavi-Jingmen Aquatic Flavivirus
# FJFL = Flavi-Jingmen Flavi-Like
# FJJI = Flavi-Jingmen Jingmenvirus
# FJUN = Flavi-Jingmen Unclassified
# 
# Pesti-LGF Clade
# PLLG = Pesti-LGF Large Genome Flavivirus
# PLPV = Pesti-LGF Pestivirus
# PLUN = Pesti-LGF Unclassified
# 
# Hepaci-Pegi Clade
# HPPV = Hepaci-Pegi Pegivirus
# HPHV = Hepaci-Pegi Hepacivirus
# HPUN = Hepaci-Pegi Unclassified
# 
# TOMB = Tombusvirus out group

clade_colors <- c(
  "Large genome Flavi-"="#DEC579",
  "Pesti-like"="#78777B",
  "Jingmen"="#82CFF1",
  "Orthoflavi-like"="#CACACA",
  "Orthoflavi"="#217B3A",
  "Hepaci"="#8C1D57",
  "Pegi"="#CC6779",
  "Pesti"="#44ABA1"
)

clade_groups <- c(
  "PLPV"="Pesti",
  "HPPV"="Pegi",
  "HPHV"="Hepaci",
  "FJJI"="Jingmen",
  "FJMB"="Orthoflavi",
  "FJNV"="Orthoflavi",
  "FJAF"="Orthoflavi",
  "FJTB"="Orthoflavi",
  "FJIS"="Orthoflavi-like",
  "FJFL"="Orthoflavi-like",
  "FJUN"="Orthoflavi-like",
  "PLLG"="Large genome Flavi-",
  "PLUN"="Pesti",
  "TOMB"="Tombusvirus"
)

# Mifsud et al. trees
m_e  <- ape::read.tree(paste(BASEDIR, "trees/refolded_fullglyco_E_3di_famsa_trim35.fas.treefile", sep=""))
m_e1 <- ape::read.tree(paste(BASEDIR, "trees/refolded_fullglyco_E1_3di_AA_famsa_parts.nex.treefile", sep=""))
m_e2 <- ape::read.tree(paste(BASEDIR, "trees/refolded_fullglyco_E2_3di_AA_famsa_parts.nex.treefile", sep=""))

# Foldmason trees
f_e_concat <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E_parts.nex.treefile", sep=""))
f_e_aa     <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E_aa_trim.fa.treefile", sep=""))
f_e_3di    <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E_3di_trim.fa.treefile", sep=""))

f_e1_concat <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E1_parts.nex.treefile", sep=""))
f_e1_aa     <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E1_aa_trim.fa.treefile", sep=""))
f_e1_3di    <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E1_3di_trim.fa.treefile", sep=""))

f_e2_concat <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E2_parts.nex.treefile", sep=""))
f_e2_aa     <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E2_aa_trim.fa.treefile", sep=""))
f_e2_3di    <- ape::read.tree(paste(BASEDIR, "trees/foldmason_E2_3di_trim.fa.treefile", sep=""))

# Compute score metrics between two trees
# Split symmetric difference metrics = Robinson-Foulds distance
compute_scores <- function(treeA, treeB) {
  quartets <- QuartetStatus(treeA, treeB)
  splits   <- SplitStatus(treeA, treeB)
  quartet_scores <- SimilarityMetrics(quartets, similarity = TRUE)
  split_scores   <- SimilarityMetrics(splits,   similarity = FALSE)
  combined <- rbind(quartet_scores, split_scores)
  row.names(combined) <- c("quartet", "split")
  return(combined) 
}

# Compute scores for trees from concat, AA-only and SS-only MSAs vs Mifsud
compute_all_scores <- function(ref, concat, aa, ss) {
  scores_concat <- compute_scores(ref, concat)
  scores_aa     <- compute_scores(ref, aa)
  scores_3di    <- compute_scores(ref, ss)
  result <- data.frame(
    tree_name = c("concat", "aa", "3di"),
    quartet_similarity = c(
      scores_concat['quartet', 'SymmetricDifference'],
      scores_aa['quartet', 'SymmetricDifference'],
      scores_3di['quartet', 'SymmetricDifference']
    ),
    rf_distance = c(
      scores_concat['split', 'SymmetricDifference'],
      scores_aa['split', 'SymmetricDifference'],
      scores_3di['split', 'SymmetricDifference']
    )
  ) 
  return(result)
}
all_scores <- bind_rows(
  compute_all_scores(m_e, f_e_concat, f_e_aa, f_e_3di)     %>% mutate(set = "E"),
  compute_all_scores(m_e1, f_e1_concat, f_e1_aa, f_e1_3di) %>% mutate(set = "E1"),
  compute_all_scores(m_e2, f_e2_concat, f_e2_aa, f_e2_3di) %>% mutate(set = "E2")
) 
all_scores_pivot <- all_scores %>%
  pivot_wider(
    names_from=tree_name,
    values_from=c(quartet_similarity, rf_distance)
  )
write.table(all_scores_pivot, paste(BASEDIR, "glycoproteins/scores.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

# Format quartet similarity and RF distance label 
get_score_label <- function(all_scores, set_, tree_name_) {
  pattern <- 'QS=%.2f\nRF=%.2f'
  row <- all_scores %>% filter(set == set_ & tree_name == tree_name_)
  text <- sprintf(pattern, row$quartet_similarity, row$rf_distance)
  return(text)
}

# Need plot limits for reliable placement of tree scale
get_plot_limits <- function(tree) {
  plot_data <- tree$data
  list(xmin = min(plot_data$x, na.rm = TRUE),
       xmax = max(plot_data$x, na.rm = TRUE),
       ymin = min(plot_data$y, na.rm = TRUE),
       ymax = max(plot_data$y, na.rm = TRUE))
}

# Main tree plotting function
plot_tree <- function(tree, all_scores, layout="ape") {
  label_clades <- as.data.frame(tree$tip.label, nm="label") %>%
    mutate(clade=clade_groups[sub("_.*$", "", label)])
  p1 <- ggtree(tree, layout=layout, size=0.2) 
  limits <- get_plot_limits(p1)
  p1 <- p1 %<+% label_clades +
    geom_tippoint(aes(fill=clade), shape=21, size=1.0, stroke=0.2, color="black", show.legend = TRUE) +
    geom_treescale(
      family = "Helvetica",
      fontsize=2,
      linesize=0.2,
      x=limits$xmin,
      y=limits$ymin,
      offset=(limits$ymax-limits$ymin)*.05,
      width=.5
    ) + 
    theme_tree() +
    coord_cartesian(clip = "off") +
    theme(text = element_text(family = "Helvetica"))
  return(p1)
}

# Plot all trees. Foldmason trees are annotated with QS/RF scores
m_e_plot  <- plot_tree(m_e)
m_e1_plot <- plot_tree(m_e1)
m_e2_plot <- plot_tree(m_e2)
f_e_plot  <- plot_tree(f_e_concat)   + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E", "concat"))
f_e_aa_plot  <- plot_tree(f_e_aa)    + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E", "aa"))
f_e_ss_plot  <- plot_tree(f_e_3di)   + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E", "3di"))
f_e1_plot <- plot_tree(f_e1_concat)  + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E1", "concat"))
f_e1_aa_plot  <- plot_tree(f_e1_aa)  + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E1", "aa"))
f_e1_ss_plot  <- plot_tree(f_e1_3di) + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E1", "3di"))
f_e2_plot <- plot_tree(f_e2_concat)  + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E2", "concat"))
f_e2_aa_plot  <- plot_tree(f_e2_aa)  + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E2", "aa"))
f_e2_ss_plot  <- plot_tree(f_e2_3di) + annotate("text", hjust=1, vjust=-.5, x=Inf, y=-Inf, size=2, label=get_score_label(all_scores, "E2", "3di"))

# Row/column titles
font <- gpar(fontsize=7, fontface='bold', fontfamily='Helvetica')
e_col_label <- wrap_elements(panel = textGrob('E', gp=font),)
e1_col_label <- wrap_elements(panel = textGrob('E1', gp=font))
e2_col_label <- wrap_elements(panel = textGrob('E2', gp=font))
m_row_label <- wrap_elements(panel = textGrob(expression(bold('Mifsud ' * bolditalic('et al.') * ' 2024')), gp=font))
f_row_label <- wrap_elements(panel = textGrob('FoldMason\nAA+3Di', gp=font))
f_aa_row_label <- wrap_elements(panel = textGrob('FoldMason\nAA', gp=font))
f_ss_row_label <- wrap_elements(panel = textGrob('FoldMason\n3Di', gp=font))

# Plot patchwork
layout <- "
#ABC
DEFG
HIJK
LMNO
PQRS
"
e_col_label + e1_col_label + e2_col_label +
  m_row_label + m_e_plot + m_e1_plot + m_e2_plot +
  f_row_label + f_e_plot + f_e1_plot + f_e2_plot +
  f_aa_row_label + f_e_aa_plot + f_e1_aa_plot + f_e2_aa_plot +
  f_ss_row_label + f_e_ss_plot + f_e1_ss_plot + f_e2_ss_plot +
  plot_layout(
    design = layout,
    widths = c(.6,.8,.8,.8),
    heights = c(.1,1,1,1,1),
    guides='collect'
  ) &
  scale_fill_manual(
    values = clade_colors,
    limits = names(clade_colors),
  ) &
  guides(fill = guide_legend(title = NULL, byrow=TRUE)) &
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size=7),
    legend.margin=margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(-.2, "lines"),
  )

ggsave(paste(BASEDIR, "glycoproteins/trees.pdf", sep=""), units="mm", width=160, height=200, dpi=300)
ggsave(paste(BASEDIR, "glycoproteins/trees.svg", sep=""), device=svg, units="mm", width=160, height=200, dpi=300)
ggsave(paste(BASEDIR, "glycoproteins/trees.png", sep=""), units="mm", width=160, height=200, dpi=300, bg="white")