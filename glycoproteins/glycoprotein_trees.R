BiocManager::install("YuLab-SMU/treedataverse")

library(ape)
library(phangorn)
library(Quartet)
library(dplyr)
library(ggtree)
library(ggplot2)
library(svglite)


# Mifsud et al trees
m_e <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/mifsud_denv1_e_mapped.nw")
m_e1 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/mifsud_bvdv1_wsv_e1_mapped.nw")
m_e2 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/mifsud_bvdv1_tdav_e2_mapped.nw")
m_e2 <- drop.tip(m_e2, c("HPHV_Hepacivirus_sp."))

# FoldMason trees
#f_e <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/foldmason_denv1_e.nw")  # 3Di tree
f_e <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/foldmason_denv1_e_aa.nw")
f_e1 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/foldmason_bvdv1_wsv_e1.nw")
f_e2 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/trees/foldmason_bvdv1_tdav_e2.nw")

## Remove OKIAV332/WHYC-1 as in their manuscript
f_e <- drop.tip(f_e, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e1 <- drop.tip(f_e1, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e2 <- drop.tip(f_e2, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))


# Refolded trees v2
m_e <- ape::read.tree("/Users/gamcil/Documents/glyco/Mifsud_et_al_Repository/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E_3di_AA_famsa_parts.nex.treefile")
m_e <- ape::read.tree("/Users/gamcil/Documents/glyco/Mifsud_et_al_Repository/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E_3di_famsa_trim35.fas.treefile")
m_e1 <- ape::read.tree("/Users/gamcil/Documents/glyco/Mifsud_et_al_Repository/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E1_3di_AA_famsa_parts.nex.treefile")
m_e2 <- ape::read.tree("/Users/gamcil/Documents/glyco/Mifsud_et_al_Repository/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E2_3di_AA_famsa_parts.nex.treefile")

f_e <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E_aa_trim35.fa.treefile")
f_e1 <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E1_aa_trim35.fa.treefile")
f_e2 <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E2_aa_trim35.fa.treefile")

f_e <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E_partitions.nex.treefile")
f_e1 <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E1_partitions.nex.treefile")
f_e2 <- ape::read.tree("/Users/gamcil/Documents/glyco/mifsud_refolded_foldmason/E2_partitions.nex.treefile")



# Function to find matching nodes and extract bootstrap values
find_matching_nodes <- function(tree1, tree2) {
  # Read the trees
  # tree1 <- read.tree(tree1_path)
  # tree2 <- read.tree(tree2_path)
  
  # Ensure both trees have the same taxa
  if (!setequal(tree1$tip.label, tree2$tip.label)) {
    stop("The trees do not have the same set of taxa.")
  }
  
  # Get all the node labels (taxa) for each internal node
  tree1_node_labels <- sapply(1:tree1$Nnode, function(i) {
    sort(tree1$tip.label[Descendants(tree1, i + length(tree1$tip.label), type = "tips")[[1]]])
  })
  tree2_node_labels <- sapply(1:tree2$Nnode, function(i) {
    sort(tree2$tip.label[Descendants(tree2, i + length(tree2$tip.label), type = "tips")[[1]]])
  })
  
  # Find matching nodes
  matching_nodes <- which(sapply(tree1_node_labels, function(node1) {
    any(sapply(tree2_node_labels, function(node2) identical(node1, node2)))
  }))
  
  # Extract bootstrap values for matching nodes
  matching_nodes_info <- data.frame(
    Tree1_Node = numeric(0),
    Tree1_Bootstrap = numeric(0),
    Tree2_Node = numeric(0),
    Tree2_Bootstrap = numeric(0)
  )
  
  for (node1 in matching_nodes) {
    node1_label <- tree1_node_labels[[node1]]
    node2 <- which(sapply(tree2_node_labels, function(node2) identical(node1_label, node2)))
    
    if (length(node2) > 0) {
      node2 <- node2[1]
      
      bootstrap1 <- tree1$node.label[node1]
      bootstrap2 <- tree2$node.label[node2]
      
      matching_nodes_info <- rbind(matching_nodes_info, data.frame(
        Tree1_Node = node1,
        Tree1_Bootstrap = as.numeric(bootstrap1),
        Tree2_Node = node2,
        Tree2_Bootstrap = as.numeric(bootstrap2)
      ))
    }
  }
  
  return(matching_nodes_info)
}

calculate_mean_tree_support <- function(tree) {
  values <- sapply(1:tree$Nnode, function(i) {
    return(as.numeric(tree$node.label[i]))
  })
  return(mean(values, na.rm=T)) 
}

calculate_matching_node_stats <- function(a, b) {
  matchingNodes <- find_matching_nodes(a, b) %>%
    mutate(difference = Tree1_Bootstrap - Tree2_Bootstrap)
  stats <- data.frame(
    Mean_Tree1_Bootstrap = calculate_mean_tree_support(a),
    Mean_Tree2_Bootstrap = calculate_mean_tree_support(b),
    Mean_Matched_Tree1_Bootstrap = mean(matchingNodes$Tree1_Bootstrap, na.rm = TRUE),
    Mean_Matched_Tree2_Bootstrap = mean(matchingNodes$Tree2_Bootstrap, na.rm = TRUE),
    Mean_Difference = mean(matchingNodes$difference, na.rm = TRUE),
    Total_Nodes = Nnode(a),
    Total_Matching_Nodes = count(matchingNodes),
    Total_Difference_Positive = sum(matchingNodes$difference > 0, na.rm = TRUE),
    Total_Difference_Negative = sum(matchingNodes$difference < 0, na.rm = TRUE),
    Total_Difference_Zero = sum(matchingNodes$difference == 0, na.rm = TRUE)
  )
  return(stats)
}


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



e_node_stats <- calculate_matching_node_stats(f_e, m_e) %>% mutate(tree="E")
e1_node_stats <- calculate_matching_node_stats(f_e1, m_e1) %>% mutate(tree="E1")
e2_node_stats <- calculate_matching_node_stats(f_e2, m_e2) %>% mutate(tree="E2")
node_stats <- rbind(e_node_stats, e1_node_stats)
node_stats <- rbind(node_stats, e2_node_stats)
node_stats

mifsud_e <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/refolded_fullglyco_E_3di_famsa_trim35.fas.treefile")
mifsud_e1 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/refolded_fullglyco_E1_3di_famsa_trim35.fas.treefile")
mifsud_e2 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/refolded_fullglyco_E2_3di_famsa_trim35.fas.treefile")
clustal_e <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/E_trimal.fa.treefile")
clustal_e1 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/E1_trimal.fa.treefile")
clustal_e2 <- read.tree("/Users/gamcil/repos/foldmason-analysis/glycoproteins/seq/E2_trimal.fa.treefile")
compute_scores(mifsud_e, clustal_e)
compute_scores(mifsud_e1, clustal_e1)
compute_scores(mifsud_e2, clustal_e2)

e_node_stats <- calculate_matching_node_stats(clustal_e, mifsud_e) %>% mutate(tree="E")
e1_node_stats <- calculate_matching_node_stats(f_e1, m_e1) %>% mutate(tree="E1")
e2_node_stats <- calculate_matching_node_stats(f_e2, m_e2) %>% mutate(tree="E2")
node_stats <- rbind(e_node_stats, e1_node_stats)
node_stats <- rbind(node_stats, e2_node_stats)
node_stats



f_e <- read.tree()  # 3Di tree
# f_e <- read.tree("trees/foldmason_denv1_e_aa.nw")
f_e1 <- read.tree("trees/foldmason_bvdv1_wsv_e1.nw")
f_e2 <- read.tree("trees/foldmason_bvdv1_tdav_e2.nw")

## Remove OKIAV332/WHYC-1 as in their manuscript
f_e <- drop.tip(f_e, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e1 <- drop.tip(f_e1, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))
f_e2 <- drop.tip(f_e2, c("FJJI_Dipteran_jingmen-related_virus_isolate_OKIAV332_segment_2_orfVP1", "FJJI_Wuhan_aphid_virus_1_strain_WHYC-1_segment_2_orfVP1"))

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
