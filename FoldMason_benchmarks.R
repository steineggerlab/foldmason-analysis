# Make plots for FoldMason manuscript
# Requires TSVs in data/:
# - homstrad_scores.tsv   - family/tool/score type/score
# - homstrad_families.tsv - family/member count
# - afdb_scores.tsv       - family/tool/score type/score
# - afdb_families.tsv     - family/CATH domains/total cluster member count/domain count/20 members analysed
# - scaling_times.tsv     - subset/tool/time
# - homstrad_lengths.csv  - protein/length
# - pdb_protein_sizes.csv - bin/count

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggrepel)
library(scales)
library(svglite)
library(patchwork)
library(scico)
library(forcats)

BASEDIR="./"

name_map <- c(
  matt = "Matt",
  mtmalign = "mTM-align",
  muscle = "MUSCLE",
  mustang = "MUSTANG",
  famsa = "FAMSA2",
  famsa3di = "FAMSA 3Di",
  caretta = "Caretta-shape", 
  mafft = "MAFFT",
  "mafft_3di" = "MAFFT 3Di",
  clustalo = "Clustal Omega",
  usalign = "US-align",
  "3dcoffee" = "3DCoffee",
  foldmason = "FoldMason",
  foldmason_fast = "FoldMason Fast",
  foldmason_refine100 = "FoldMason R100",
  foldmason_refine100_rosc = "FoldMason R100 OSC",
  homstrad = "HOMSTRAD",
  "foldseek" = "Foldseek",
  "foldseektm" = "Foldseek-TM",
  "clesw" = "CLE-SW",
  "ce" = "CE",
  "dali" = "DALI",
  "tmalign" = "TM-align",
  "mmseqs" = "MMseqs2",
  "fm_probs" = "FoldMason Probs"
)

colour_map <- c(
  "FoldMason" = "#FF0000",
  "FoldMason Fast" = "pink",
  "FoldMason R100" = "#AA0099",
  "3DCoffee" = "goldenrod",
  "Caretta-shape" = "#D55E00",
  "MAFFT 3Di" = "blue",
  "MUSTANG" = "darkgreen",
  "Matt" = "#CC79A7",
  "US-align" = "black",
  "mTM-align" = "#56B4E9",
  "Clustal Omega" = "#E69F00",
  "FAMSA2" = "#0072B2",
  "MAFFT" = "#999999",
  "MUSCLE" = "#009E73"
)

shape_map <- c(
  "Sequence-based" = 1,
  "Structure-based" = 2,
  "FoldMason" = 17
)


common_theme <- theme_bw(base_family="Helvetica", base_size=7) + theme(
  panel.border = element_blank(),
  axis.line = element_line(linewidth = unit(0.2, 'pt')),
  axis.title.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  strip.background = element_blank(),
  plot.margin = margin(1, 1, 1, 1),
  legend.position="bottom",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.title=element_blank(),
  legend.key=element_blank(),
  legend.key.size=unit(0.6, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
)

structure_tools = c("mTM-align", "MUSTANG", "Matt", "Caretta-shape",
                    "FoldMason", "FoldMason R100", "FoldMason Fast", "US-align", "3DCoffee",
                    "FoldMason R100 OSC", "MAFFT 3Di")
sequence_tools = c("Clustal Omega", "MUSCLE", "MAFFT", "FAMSA2")

# Panel 1: Homstrad
homstrad = new.env()

# 55 hard families (<25% identity, >=4 members)
homstrad.hard <- c(
  "ABC_tran", "Acetyltransf", "alpha-amylase", "alpha-amylase_C", "alpha-amylase_NC", "Asp_Glu_race_D",
  "blmb", "bv", "cat3", "CH", "CPSase_L_chain", "DEATH", "DHOdehase", "Epimerase", "FAD_binding_4",
  "FAD-oxidase_C", "FAD-oxidase_NC", "fn3", "GEL", "ghf18", "ghf33", "ghf5", "Glyco_hydro_18_D2",
  "Haloperoxidase", "HATPase_c", "helicase_C", "HGTP_anticodon", "histone", "HMG_box", "hormone_rec",
  "igI", "igV", "int", "kinase", "laminin_G", "lipocalin", "oat", "p450", "pdc", "Peptidase_M24",
  "PH", "Phage_integrase", "porin", "proteasome", "reductases", "Rhodanese", "Ribosomal_L6_D", "rrm",
  "sdr", "Sm", "sugbp", "Sulfotransfer", "thiored", "TPR", "tRNA-synt_2b" 
)

homstrad.data <- read.delim(
  paste(BASEDIR, "data/homstrad_scores.tsv", sep=""),
  sep="\t",
  header = F,
  col.names = c("family", "tool", "type", "score"),
) %>%
  mutate(
    osc=str_ends(tool, "_osc"),
    tool=name_map[if_else(osc, str_remove(tool, "_osc$"), tool)],
    base=if_else(tool %in% structure_tools, "Structure-based", "Sequence-based")
  )

homstrad.counts <- read.delim(
  paste(BASEDIR, "data/homstrad_families.tsv", sep=""),
  sep="\t",
  header=F,
  col.names=c("family", "count")
)
 
# Sum-of-pairs plot
# Pairs from reference MSA in test MSA (specificity)
# Pairs from test MSA in reference MSA (accuracy)
homstrad.sop_scores <- homstrad.data %>%
  filter(tool != "FoldMason R100") %>%
  left_join(homstrad.counts, by=join_by(family)) %>%
  filter(count >= 2 & startsWith(type, "sp_") & tool != "Homstrad") %>%
  mutate(score=as.numeric(score), tool=if_else(tool == "FoldMason R100 OSC", "FoldMason R100", tool)) %>%
  group_by(tool, type, base) %>%
  summarise(mean=mean(score), sd=sd(score), .groups="drop") %>%
  pivot_wider(names_from=type, values_from=c(mean, sd), names_sep="_") %>%
  mutate(f1 = 2 * (mean_sp_fwd * mean_sp_rev) / (mean_sp_fwd + mean_sp_rev), subset="HOMSTRAD Full")

homstrad.sop_scores_hard <- homstrad.data %>%
  filter(family %in% homstrad.hard & tool != "FoldMason R100") %>%
  left_join(homstrad.counts, by=join_by(family)) %>%
  filter(count >= 2 & startsWith(type, "sp_") & tool != "Homstrad") %>%
  mutate(score=as.numeric(score), tool=if_else(tool == "FoldMason R100 OSC", "FoldMason R100", tool)) %>%
  group_by(tool, type, base) %>%
  summarise(mean=mean(score), sd=sd(score), .groups="drop") %>%
  pivot_wider(names_from=type, values_from=c(mean, sd), names_sep="_") %>%
  mutate(f1 = 2 * (mean_sp_fwd * mean_sp_rev) / (mean_sp_fwd + mean_sp_rev), subset="HOMSTRAD Hard")


# Display order and shape scale, show FoldMason points filled and other tools hollow
homstrad.combined <- bind_rows(homstrad.sop_scores, homstrad.sop_scores_hard)
homstrad.combined$plot_shape <- as.character(homstrad.combined$base)
homstrad.combined <- homstrad.combined %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
homstrad.sop_scores$plot_shape <- as.character(homstrad.sop_scores$base)
homstrad.sop_scores <- homstrad.sop_scores %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
tool_levels <- names(colour_map)
homstrad.combined$tool <- factor(homstrad.combined$tool, levels = tool_levels)
level_to_shape <- unique(homstrad.combined[, c("tool", "plot_shape")])
level_to_shape <- level_to_shape[match(tool_levels, level_to_shape$tool), ]
correct_shapes <- unname(shape_map[level_to_shape$plot_shape])

breaks_fun <- function(x) {
  if (min(x) < 65) {
    seq(55, 100, 10)
  } else {
    seq(70, 100, 5)
  }
}
breaks_fun_y <- function(x) {
  if (min(x) < 65) {
    seq(50, 100, 10)
  } else {
    seq(80, 100, 5)
  }
}

# Fig. 2A
homstrad.plot <- ggplot(homstrad.combined) +
  aes(y=mean_sp_fwd, x=mean_sp_rev, colour=tool, shape=plot_shape) +
  facet_wrap(~ subset, nrow=2, scales="free") +
  scale_color_manual(values = colour_map) +
  scale_shape_manual(values = shape_map, breaks = c("Structure-based", "Sequence-based")) +
  geom_point(size=0.8, stroke=0.3) +
  scale_y_continuous(breaks = breaks_fun_y) +
  scale_x_continuous(breaks = breaks_fun) +
  guides(
    shape=guide_legend(nrow=2),
    color=guide_legend(
      nrow=2,
      override.aes = list(shape = correct_shapes, linetype = 0)
    ),
  ) +
  labs(x="Sensitivity (%)", y="Specificity (%)", colour="Tool", shape="Input type") +
  common_theme

homstrad.plot
ggsave(file=paste(BASEDIR, "figures/Fig_2A_HOMSTRAD_SoP.png", sep=""), units="mm", width=140, height=140, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_2A_HOMSTRAD_SoP.pdf", sep=""), units="mm", width=140, height=140, dpi=300, bg="white")

# F1 scores
homstrad.sop_scores %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1) %>%
  arrange(desc(f1))
homstrad.sop_scores_hard %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1) %>%
  arrange(desc(f1))

# Structure-based vs sequence-based tools
homstrad.sop_scores %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1) %>%
  group_by(base) %>%
  summarise(mean_sp_fwd=mean(mean_sp_fwd), mean_sp_rev=mean(mean_sp_rev), mean_f1=mean(f1))
homstrad.sop_scores_hard %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1) %>%
  group_by(base) %>%
  summarise(mean_sp_fwd=mean(mean_sp_fwd), mean_sp_rev=mean(mean_sp_rev), mean_f1=mean(f1))

# All metric boxplots
homstrad.other_metrics <- homstrad.data %>%
  filter(!osc & type %in% c("lddt", "sp_fwd", "sp_rev")) %>%
  left_join(homstrad.counts, by=join_by(family)) %>% 
  filter(count >= 4) %>%
  mutate(score=as.numeric(score))
homstrad.other_metrics$tool <- factor(homstrad.other_metrics$tool, levels=unique(homstrad.other_metrics$tool[order(homstrad.other_metrics$base, homstrad.other_metrics$tool)]))
ggplot(homstrad.other_metrics) +
  aes(y=score, x=tool, fill=base) +
  facet_wrap(~type, scales = "free") +
  geom_boxplot() +
  common_theme +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
ggsave(file=paste(BASEDIR, "figures/homstrad_all_metrics.pdf", sep=""), units="mm", width=300, height=300, dpi=300, bg="white")

# Fig. S3
# Correlation between LDDT and SoP reference-based scores on HOMSTRAD
homstrad.wide <- homstrad.data %>%
  filter(!osc & tool != "HOMSTRAD" & type %in% c("lddt", "sp_fwd", "sp_rev")) %>%
  pivot_wider(names_from=type, values_from=score, values_fn=as.numeric) %>%
  left_join(homstrad.counts, by=c("family")) %>%
  filter(count >= 4)
homstrad.longer <- homstrad.wide %>%
  pivot_longer(cols=c(sp_fwd, sp_rev), names_to="score_type", values_to="score") %>%
  mutate(score_type = factor(score_type, levels=c('sp_fwd', 'sp_rev')))
homstrad.correlations <- homstrad.longer %>% 
  filter(count >= 4) %>%
  group_by(score_type) %>%
  summarize(
    pearson=cor(x=score, y=lddt, method="pearson"),
    spearman=cor(x=score, y=lddt, method="spearman")
  )
ggplot(homstrad.longer %>% filter(count >= 4)) +
  aes(x=lddt, y=score, colour=score_type) + 
  facet_wrap(
    ~score_type,
    nrow=1,
    strip.position = "left",
    labeller = as_labeller(c(
      sp_fwd = "Sum-of-pairs forward (Sensitivity)",
      sp_rev = "Sum-of-pairs reverse (Specificity)"
    )),
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth=0.5, color="black") +
  geom_point(size=0.4, alpha=.5) +
  geom_text(
    data = homstrad.correlations,
    aes(x=Inf, y=50, label = paste( "Spearman R =", round(spearman, 2), "\nPearson R =", round(pearson, 2))),
    color="black",
    hjust = 1.1,
    vjust = 5.0,
    size = 2,
  ) +
  labs(
    x="Average MSTA LDDT",
    y=NULL
  ) +
  guides(shape="none", colour="none",) +
  coord_fixed(ratio=0.01) +
  common_theme +
  theme(
    strip.background=element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=6)
  )
ggsave(file=paste(BASEDIR, "figures/Fig_S3_HOMSTRAD_correlations.pdf", sep=""), units="mm", width=140, height=80, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_S3_HOMSTRAD_correlations.png", sep=""), units="mm", width=140, height=80, dpi=300, bg="white")


# Panel 2: 1000 AFDB Clusters
afdb = new.env()

# 40 clusters used for parameter optimisation, removed for benchmark
afdb.training <- c(
  "S7TGR0", "R9KI14", "M5RUD2", "L9WUI5", "K6ACG0", "K0EYU6", "I4F5N6", "E8N900",
  "D3TB77", "B4WMN9", "W1PFB2", "R7ETQ2", "R1FSM3", "L8LTQ4", "K2G3Q9", "F4NU11",
  "C3XV33", "B9LR42", "B5GRN6", "B2A5F0", "W7QDU7", "W0E3L2", "V2YI92", "U7L1T7",
  "T1FNF4", "R8GWM8", "R8APS3", "R6WAH6", "Q3SL35", "Q384X5", "A0A838SJH4",
  "A0A841HGH3", "A0A842R2G4", "A0A846H7L1", "A0A847S9B0", "A0A852XNB6", "B2IJL3",
  "C7LW54", "C9LSN8", "F6CXK8" 
)
afdb.data_base <- read.delim(
  paste(BASEDIR, "data/afdb_scores.tsv", sep=""),
  sep="\t",
  header=FALSE,
  col.names=c("family", "tool", "type", "score")
) %>%
  filter(type != "lddt_columns" & !(family %in% afdb.training)) %>%
  mutate(
    osc=str_ends(tool, "_osc"),
    tool=if_else(osc, str_remove(tool, "_osc$"), tool),
    base=ifelse(name_map[tool] %in% structure_tools, "Structure-based", "Sequence-based"),
    score=as.numeric(score)
  )

# Cluster statistics, need num_domains
afdb.family_data <- read.delim(
  paste(BASEDIR, "data/afdb_families.tsv", sep=""),
  sep="\t",
  header=FALSE,
  col.names=c("rep_accession", "domains", "num_members", "num_domains", "member_accessions")
) %>%
  mutate(num_domains = as.factor(num_domains))

afdb.data_base <- afdb.data_base %>%
  left_join(afdb.family_data, by=c('family' = 'rep_accession'))

# Structure vs sequence-based tools
# LDDT scores
afdb.data_base %>%
  filter(type == "lddt" & osc) %>%
  group_by(base) %>%
  summarise(avg_lddt=mean(score, na.rm=TRUE)) %>%
  pivot_wider(names_from = base, values_from = avg_lddt) %>%
  mutate(difference = (`Structure-based` - `Sequence-based`) * 100)

# Runtimes
afdb.data_base %>%
  filter(type == "time") %>%
  group_by(base) %>%
  summarise(avg_lddt=mean(score, na.rm=TRUE)) %>%
  pivot_wider(names_from = base, values_from = avg_lddt) %>%
  mutate(difference = (`Structure-based` / `Sequence-based`))

# LDDT improvement by refinement
afdb.refine_diff <- afdb.data_base %>%
  filter(tool %in% c("foldmason_fast", "foldmason_refine100") & type == "lddt" & osc) %>%
  dplyr::select(family, tool, score) %>%
  pivot_wider(names_from=tool, values_from=score) %>%
  mutate(delta = `foldmason_refine100` - `foldmason_fast`)
afdb.refine_diff %>%
  summarise(mean=mean(delta), median=median(delta), min_diff=min(delta) * 100, max_diff=max(delta) * 100)

# Calculate tool speedup
afdb.summary_time <- afdb.data_base %>%
  filter(type == "time") %>%
  group_by(tool) %>%
  summarise(
    avg_time=mean(score),
    min_time=min(score),
    max_time=max(score)
  ) %>%
  mutate(
    speedup=max(avg_time) / avg_time,
    slowdown=avg_time / min(avg_time),
  )

# Summarise full-MSA LDDT scores
summariseData <- function(scoreDf, timeDf) {
  summary <- scoreDf %>%
    filter(type == "lddt") %>%
    group_by(tool) %>%
    summarise(
      avg_lddt=mean(score),
      avg_lddt_1dom=mean(score[num_domains==1]),
      avg_lddt_2dom=mean(score[num_domains==2]),
      avg_lddt_3dom=mean(score[num_domains==3]),
      avg_lddt_4dom=mean(score[num_domains==4]),
      min_score=min(score),
      max_score=max(score),
      percentile_10 = quantile(score, 0.1, na.rm = TRUE),
      percentile_90 = quantile(score, 0.90),
      first_quartile = quantile(score, 0.25),
      third_quartile = quantile(score, 0.75),
      stdev=sd(score),
    ) %>%
    left_join(timeDf, by = c("tool")) %>%
    mutate(
      tool=name_map[tool],
      base=if_else(tool %in% structure_tools, "Structure-based", "Sequence-based"),
    )
  fm_only      <- summary %>% dplyr::filter(tool == "FoldMason")
  fm_r100_only <- summary %>% dplyr::filter(tool == "FoldMason R100")
  return(
    summary %>%
      mutate(
        fm_lddt_diff = (avg_lddt - fm_only %>% pull(avg_lddt)) * 100,  # difference to foldmason
        fm_r100_lddt_diff = (avg_lddt - fm_r100_only %>% pull(avg_lddt)) * 100,  # difference to refined foldmason
        fm_speedup = avg_time / fm_only %>% pull(avg_time),   # speedup wrt foldmason
        fm_slowdown = fm_only %>% pull(avg_time) / avg_time,  # slowdown wrt foldmason
      )
  )
}
afdb.summary <- summariseData(afdb.data_base %>% filter(!osc), afdb.summary_time) %>% filter(!str_starts(tool, "FoldMason R100"))
afdb.summary %>%
  dplyr::select(tool, avg_lddt, avg_lddt_1dom, avg_lddt_2dom, avg_lddt_3dom, avg_lddt_4dom, avg_time)
afdb.summary_osc <- summariseData(afdb.data_base %>% filter(osc), afdb.summary_time) %>% filter(!str_starts(tool, "FoldMason R100"))
afdb.summary_osc %>%
  dplyr::select(tool, avg_lddt, avg_lddt_1dom, avg_lddt_2dom, avg_lddt_3dom, avg_lddt_4dom, avg_time)

# Fig. S7
# Per-tool performance boxplots on AFDB
ggplot(
  afdb.data_base %>%
    filter(tool != "foldmason_refine100" & tool != "foldmason_refine100_rosc" & type == "lddt") %>%
    mutate(tool=name_map[tool])
  ) +
  aes(x=reorder(tool, score, FUN=median), color=tool, y=score) +
  facet_wrap(
    ~ osc,
    nrow=2,
    labeller = labeller(osc = c("TRUE" = "Only scoring columns normalized LDDT", "FALSE" = "Average MSTA LDDT"))) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.2), alpha=0.2) +
  labs(x="Tool", y="Average MSTA LDDT") +
  scale_color_manual(values = colour_map) +
  scale_x_discrete(guide=guide_axis(angle=45)) +
  guides(color="none") +
  common_theme
ggsave(file=paste(BASEDIR, "figures/Fig_S7_AFDB_all_tool_boxplots.png", sep=""), units="mm", limitsize=F, width=140, height=120, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_S7_AFDB_all_tool_boxplots.pdf", sep=""), units="mm", limitsize=F, width=140, height=120, bg="white")


# Fig. 2B
# Average LDDT over AFDB
afdb.summary$plot_shape <- as.character(afdb.summary$base)
afdb.summary <- afdb.summary %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
afdb.summary_osc$plot_shape <- as.character(afdb.summary_osc$base)
afdb.summary_osc <- afdb.summary_osc %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
line_test_data <- afdb.data_base %>%
    filter(type == 'lddt' & !osc & !str_starts(tool, "foldmason_refine100")) %>%
    mutate(tool=name_map[tool], base=if_else(tool %in% structure_tools, "Structure-based", "Sequence-based"))
line_samples <- line_test_data %>%
  group_by(num_domains, family) %>%
  summarize(med = median(score), .groups='drop') %>%
  arrange(num_domains, med) %>%
  group_by(num_domains) %>%
  mutate(score = med, row = row_number()) %>%
  filter(row %% 10 == 1)  # Sort by median, take every 10th structure

# Get avg score per domain count level
overall_avg_data <- line_test_data %>%
  group_by(num_domains, tool) %>%
  summarize(score = mean(score), .groups='drop') %>%
  mutate(family = "Average", base = if_else(tool %in% structure_tools, "Structure-based", "Sequence-based"))

line_test_data_sampled <- line_test_data %>%
  semi_join(line_samples, by=c("num_domains", "family"))

line_test_data_combined <- bind_rows(line_test_data_sampled, overall_avg_data)

line_test_levels_sampled <- line_test_data_sampled %>%
  group_by(family) %>%
  summarize(med = median(score), .groups='drop') %>%
  arrange(med) %>%
  pull(family)

line_test_levels <- c("Average", line_test_levels_sampled)

line_test_data_combined <- line_test_data_combined %>%
  mutate(family = factor(family, levels=line_test_levels))

# Only show reduced set in main plot; remove filter when generating Figs. S5-6 
line_test_subset <- line_test_data_combined  %>%
  filter(tool %in% c("FoldMason", "FoldMason Fast", "3DCoffee", "FAMSA2", "MAFFT", "MAFFT 3Di", "US-align", "MUSTANG"))
line_test_data_combined <- line_test_subset %>%
  mutate(nudge_x=case_when(str_starts(tool, "FoldMason") ~ 0.0, base == "Sequence-based" ~ -0.15, base == "Structure-based" ~ 0.15))


line_test_subset$plot_shape <- as.character(line_test_subset$base)
line_test_subset <- line_test_subset %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))

afdb.line_plot <- ggplot(line_test_subset) +
  aes(x=family, y=score, color=tool, shape=plot_shape) +
  facet_grid(cols=vars(num_domains), scales="free_x", labeller = labeller(num_domains = ~ paste(.x, "domain"))) +
  geom_line(aes(group = tool), linetype="dotted", linewidth=0.1, show.legend = F) +
  geom_point(size=1, stroke = .2, position=position_nudge(x=line_test_data_combined$nudge_x)) + #, position=position_dodge(width=0.5, preserve="single", )) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_shape_manual(values = shape_map, breaks = c("Structure-based", "Sequence-based")) +
  scale_color_manual(values = colour_map) +
  guides(shape="none", color="none") +
  # guides(
  #   # x="none",
  #   shape=guide_legend(nrow=4),
  #   # shape="none",
  #   color=guide_legend(
  #     nrow=4,
  #     override.aes = list(shape = shape_map[afdb.summary$plot_shape], linetype = 0)
  #   ),
  # ) +
  labs(x="Multi-domain AFDB Protein family", y="Average MSTA LDDT") +
  common_theme +
  theme(
    panel.grid.major.x = element_line(color="gray90", linewidth=0.2),
    axis.text.x = element_text(angle=45, hjust=1, size = unit(4, "pt")),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size=5,),
    axis.title.y = element_text(size=5),
    axis.text.y = element_text(size=4),
    strip.text = element_text(size=5),
    legend.text = element_text(size=5),
    legend.key.spacing.y = unit(0.01, "cm"),
    legend.key.spacing.x = unit(0.01, "cm"),
    plot.margin = margin(t=0, r=0, b=0, l=0),
    legend.box.margin = margin(t=-5, r=0, b=0, l=0),
    legend.margin = margin(t=0, r=0, b=0, l=0)
  )
afdb.line_plot
ggsave(file=paste(BASEDIR, "figures/Fig_2B_AFDB.png", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_2B_AFDB.pdf", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")

# When including all tools
ggsave(file=paste(BASEDIR, "figures/Fig_S5_AFDB_all_tools.png", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_S5_AFDB_all_tools.pdf", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")

# When using only-scoring-cols LDDT metric
ggsave(file=paste(BASEDIR, "figures/Fig_S6_AFDB_osc.png", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_S6_AFDB_osc.pdf", sep=""), units="mm", limitsize=F, width=105, height=70, bg="white")

# Breakdown of all tool LDDT scores on all AFDB clusters
family_order <- afdb.data_base %>%
  filter(type == "lddt", tool == "foldmason") %>%
  group_by(family) %>%
  summarise(ordering_score = mean(score), .groups = "drop") %>%
  arrange(desc(ordering_score)) %>%
  pull(family)
afdb.all_fams <- afdb.data_base %>% filter(type == "lddt") %>% group_by(family, tool, num_domains, base) %>% summarise(avg=mean(score), .groups='drop') %>% group_by(num_domains) %>% mutate(tool=name_map[tool], family = factor(family, levels = family_order))
ggplot(afdb.data_base %>% filter(type == "lddt") %>% group_by(family, tool, num_domains, base) %>% summarise(avg=mean(score), .groups='drop') %>% group_by(num_domains) %>% mutate(tool=name_map[tool], family = factor(family, levels = family_order))) +
  facet_wrap(vars(num_domains), dir='v', ncol=1, scales='free_y') +
  scale_y_discrete(drop = TRUE) +
  aes(y=family, x=avg, color=tool, group=tool, shape=base) +
  scale_color_manual(name="tool", values = colour_map) +
  geom_point(size=3) +
  theme_minimal()
ggsave(file=paste(BASEDIR, "figures/afdb_all_families.pdf", sep=""), units="mm", limitsize=F, width=300, height=1000, bg="white")


# Panel 3: Speed benchmark
speed = new.env()
speed.times <- read.delim(
  paste(BASEDIR, "data/afdb_scaling_times.tsv", sep=""),
  sep="\t",
  header=FALSE,
  col.names=c("subset", "tool", "type", "time")
  ) %>%
  mutate(tool=name_map[tool], base=if_else(tool %in% structure_tools, "Structure-based", "Sequence-based"))
speed.times$subsetLabel <- gsub("subset_", "", speed.times$subset)
speed.times$subsetNum <- as.numeric(speed.times$subsetLabel)
speed.times <- speed.times[order(speed.times$subsetNum),]
speed.times$subsetFactor <- factor(speed.times$subsetLabel, levels=unique(speed.times$subsetLabel))
speed.times %>%
  dplyr::select(subset, tool, time) %>%
  pivot_wider(names_from=tool, values_from=time) %>%
  dplyr::select(subset, `3DCoffee`, FAMSA2, FoldMason, `FoldMason Fast`)

convert_time <- function(seconds) {
  sapply(seconds, function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x < 60) {
      return(paste(x, "sec"))
    } else if (x < 3600) {
      return(paste(round(x / 60, 1), "min"))
    } else if (x < 86399) {
      return(paste(round(x / 3600, 1), "hr"))
    } else {
      return(paste(round(x / 86400, 1), "day"))
    }
  })
}
scientific_10 <- trans_format("log10", math_format(10^.x))
speed.times$plot_shape <- as.character(speed.times$base)
speed.times <- speed.times %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
options(scipen=999)

# Fig. 2C
# Scaling benchmark on large AFDB cluster
speed_plot <- ggplot(speed.times %>% filter(subsetNum %in% c(10, 100, 1000, 10000, 100000))) +
  aes(x=subsetNum, y=time, color=tool, group=tool, shape=plot_shape) +
  geom_point(size=1.2) +
  geom_line(alpha=0.5, linetype='dotted') +
  scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000), labels=scientific_10) +
  scale_y_log10(breaks=c(1, 10, 60, 600, 3600, 18000, 86400), labels=convert_time) +
  scale_color_manual(name="tool", values = colour_map) +
  scale_shape_manual(name="tool", values = shape_map) +
  labs(x="AFDB Subset size", y="Time (s)", color="Tool", tag="iii)") +
  guides(shape="none", colour="none") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
speed_plot
ggsave(file=paste(BASEDIR, "figures/Fig_2C_AFDB_scaling.png", sep=""), units="mm", width=80, height=80, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_2C_AFDB_scaling.pdf", sep=""), units="mm", width=80, height=80, dpi=300, bg="white")


# Fig. 2D
# Performance on flexible proteins
flex = new.env()
flex.data <- read.delim(
  paste(BASEDIR, "data/flexibility_scores.tsv", sep=""),
  sep="\t",
  header=FALSE,
  col.names=c("protein", "tool", "type", "lddt")
) %>%
  mutate(
    tool=name_map[tool],
    base=if_else(tool %in% structure_tools, "Structure-based", "Sequence-based")
  )
flex.data$plot_shape <- as.character(flex.data$base)
flex.data <- flex.data %>% mutate(plot_shape = if_else(str_starts(tool, "FoldMason"), "FoldMason", base))
flex.plot <- ggplot(flex.data) +
  aes(x=protein, y=lddt, color=tool, shape=plot_shape) +
  scale_color_manual(name="tool", values = colour_map) +
  scale_shape_manual(values = shape_map) +
  scale_x_discrete(labels=c("Calmodulin", "Globin", "OBP")) +
  ylim(0, 1) +
  facet_wrap(~ type, nrow=2, scales="free_y", labeller=labeller(type = c("msa" = "MSTA", "pair" = "Pairwise"))) + 
  guides(color="none", shape="none") +
  geom_point(position=position_jitter(width=0.3), stat="identity") +
  labs(x = "Flexible Protein", y = "Average MSTA LDDT") +
  # geom_bar(position=position_dodge(width=0.8), stat="identity", color="black", linewidth=0.2) +
  common_theme +
  theme(
    panel.grid.major.x = element_line(color="gray90", linewidth=0.2),
    axis.text.x = element_text(angle=45, hjust=1, size = unit(4, "pt")),
    axis.ticks.x = element_blank(),
    strip.placement = "inside",
  )

flex.plot
ggsave(file=paste(BASEDIR, "figures/Fig_2D_flexible_proteins.png", sep=""), units="mm", width=80, height=80, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_2D_flexible_proteins.pdf", sep=""), units="mm", width=80, height=80, dpi=300, bg="white")


# Plot entire Fig. 2
(homstrad.plot + afdb.line_plot + speed_plot + flex.plot) +
  plot_layout(axes="collect", guides="collect", widths=c(1, 3, 1, 1)) +
  plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=theme(plot.tag=)) &
  theme(
    legend.position="bottom",
    plot.tag=element_text(face="bold"),
    axis.ticks = element_line(linewidth = unit(0.2, 'pt')),
    axis.title.x = element_text(size=5,),
    axis.title.y = element_text(size=5),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4),
    strip.text = element_text(size=5),
    legend.text = element_text(size=5),
    legend.key.spacing.y = unit(0.01, "cm"),
    legend.key.spacing.x = unit(0.01, "cm"),
    plot.margin = margin(t=0, r=0, b=0, l=0),
    legend.box.margin = margin(t=-5, r=0, b=0, l=0),
    legend.margin = margin(t=0, r=0, b=0, l=0)
  )

# device=svg is required for correct import into figma, converts text to paths though
# svglite preserves text but has weirdness with path stroke/fill when imported into figma
ggsave(file=paste(BASEDIR, "figures/Fig_2_full.pdf", sep=""), units="mm", width=160, height=60, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_2_full.svg", sep=""), device=svg, units="mm", width=160, height=60, dpi=300, bg="white")


# Fig. S2
# PDB vs Homstrad protein size distributions
pdb_protein_size <- read.delim(
  paste(BASEDIR, "data/pdb_protein_sizes.csv", sep=""),
  sep=",",
  col.names = c("range", "count"),
  header = FALSE
) %>% mutate(frequency = count / sum(count), source="PDB")
homstrad_protein_size <- read.delim(
  paste(BASEDIR, "data/homstrad_lengths.csv", sep=""),
  sep=",",
  col.names = c("name", "length")
)
bins <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, Inf)
homstrad_protein_size$range <- cut(
  homstrad_protein_size$length,
  breaks=bins,
  labels=pdb_protein_size$range,
  right=FALSE
)
homstrad_protein_ranges <- as.data.frame(table(homstrad_protein_size$range))
colnames(homstrad_protein_ranges) <- c("range", "count")
homstrad_protein_ranges <- homstrad_protein_ranges %>%
  mutate(frequency = count / sum(count), source="Homstrad")
mean(homstrad_protein_size$length)
median(homstrad_protein_size$length)
combined <- rbind(pdb_protein_size, homstrad_protein_ranges)
combined$range <- factor(combined$range, levels=pdb_protein_size$range, ordered=TRUE)
ggplot(combined, aes(x = range, y=frequency, fill=source)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Protein Length (AA)", y = "Frequency", fill = "Source") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(0, 0, 0, 0),
    legend.position=c(0.9, 0.9)
  )
ggsave(file=paste(BASEDIR, "figures/Fig_S2_HOMSTRAD_PDB_distribution.png", sep=""), units="mm", width=100, height=80, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/Fig_S2_HOMSTRAD_PDB_distribution.pdf", sep=""), units="mm", width=100, height=80, dpi=300, bg="white")


# Total proteins in each database
sum(pdb_protein_size$count)
sum(homstrad_protein_ranges$count)

# AFDB per-tool MSA length distributions
afdb.lengths <- read.csv(
  paste(BASEDIR, "data/afdb_lengths.tsv", sep = ""),
  sep = "\t",
  col.names=c("family", "tool", "length")
) %>%
  mutate(tool=name_map[tool]) %>%
  left_join(afdb.family_data %>% dplyr::select(rep_accession, num_domains), by=c('family' = 'rep_accession'))

afdb.lengths_all <- afdb.lengths %>% mutate(num_domains = "All")
afdb.lengths <- rbind(afdb.lengths, afdb.lengths_all)
afdb.lengths$tool <- as.factor(afdb.lengths$tool)

ggplot(afdb.lengths) +
  aes(x = reorder(tool, length, FUN=median), y = length, color=tool) +
  scale_color_manual(name="tool", values = colour_map) +
  geom_boxplot() +
  facet_grid(rows = vars(num_domains), scales="free_y", space="fixed") +
  common_theme
ggsave(file=paste(BASEDIR, "figures/afdb_length_boxplots.pdf", sep=""), units="mm", width=180, height=160, dpi=300, bg="white")