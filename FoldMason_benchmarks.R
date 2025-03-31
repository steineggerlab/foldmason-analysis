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

BASEDIR="./foldmason-analysis/"

name_map <- c(
  matt = "Matt",
  mtmalign = "mTM-align",
  muscle = "MUSCLE",
  mustang = "MUSTANG",
  famsa = "FAMSA",
  famsa3di = "FAMSA 3Di",
  caretta = "Caretta", 
  mafft = "MAFFT",
  clustalo = "Clustal Omega",
  usalign = "US-align",
  "3dcoffee" = "3D-Coffee",
  foldmason = "FoldMason",
  foldmason_refine100 = "FoldMason R100",
  homstrad = "HOMSTRAD",
  "foldseek" = "Foldseek",
  "foldseektm" = "Foldseek-TM",
  "clesw" = "CLE-SW",
  "ce" = "CE",
  "dali" = "DALI",
  "tmalign" = "TM-align",
  "mmseqs" = "MMseqs2"
)

colour_map <- c(
  "Clustal Omega" = "#E69F00",
  "mTM-align" = "#56B4E9",
  "US-align" = "black",
  "MUSCLE" = "#009E73",
  "MUSTANG" = "#B15928",
  "FAMSA" = "#0072B2",
  "Caretta" = "#D55E00",
  "Matt" = "#CC79A7",
  "MAFFT" = "#999999",
  "3D-Coffee" = "gold",
  "FoldMason" = "#FF0000",
  "FoldMason R100" = "#AA0099"
)

shape_map <- c(
  "Sequence-based" = 16,
  "Structure-based" = 17
)

common_theme <- theme_bw(base_family="Helvetica", base_size=7) + theme(
  plot.margin = margin(1, 1, 1, 1),
  legend.position="bottom",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.title=element_blank(),
  legend.key=element_blank(),
  legend.key.size=unit(0.6, "lines"),
  legend.margin = margin(0, 0, 0, 0), # Minimize margin around the legend
  legend.box.margin = margin(0, 0, 0, 0), # Minimize box margin  
)

structure_tools = c("mTM-align", "MUSTANG", "Matt", "Caretta",
                    "FoldMason", "FoldMason New", "FoldMason New R1000", "FoldMason FWBW", "FoldMason FWBW R1000", "FoldMason R100", "US-align", "3D-Coffee")
sequence_tools = c("Clustal Omega", "MUSCLE", "MAFFT", "FAMSA")

# Panel 1: Homstrad
homstrad = new.env()
homstrad.data <- read.delim(
  paste(BASEDIR, "data/homstrad_scores.tsv", sep=""),
  sep="\t",
  header = F,
  col.names = c("family", "tool", "type", "score"),
) %>%
  mutate(
    tool=name_map[tool],
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
  left_join(homstrad.counts, by=join_by(family)) %>% 
  filter(count >= 2 & startsWith(type, "sp_") & tool != "Homstrad") %>%
  mutate(score=as.numeric(score)) %>%
  group_by(tool, type, base) %>%
  summarise(mean=mean(score), sd=sd(score), .groups="drop") %>%
  pivot_wider(names_from=type, values_from=c(mean, sd), names_sep="_") %>%
  mutate(f1 = 2 * (mean_sp_fwd * mean_sp_rev) / (mean_sp_fwd + mean_sp_rev))

homstrad.sop_scores %>% arrange(desc(f1))

homstrad.plot <- ggplot(homstrad.sop_scores) +
  aes(y=mean_sp_fwd, x=mean_sp_rev, colour=tool, shape=base) +
  xlim(78, 95) +
  ylim(78, 95) +
  scale_color_manual(values = colour_map) +
  geom_point(size=1.2) +
  guides(
    shape="none",
    color="none"
  ) +
  labs(x="Sensitivity (%)", y="Specificity (%)", colour="Tool", shape="Input type", tag="i)") +
  common_theme

homstrad.plot

ggsave(file=paste(BASEDIR, "figures/homstrad_sop_plot.pdf", sep=""), units="mm", width=100, height=100, dpi=300, bg="white")

# F1 scores
homstrad.sop_scores %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1)
homstrad.sop_scores %>%
  #dplyr::filter(tool != "Caretta") %>%
  dplyr::select(tool, base, mean_sp_fwd, mean_sp_rev, f1) %>%
  group_by(base) %>%
  summarise(mean_sp_fwd=mean(mean_sp_fwd), mean_sp_rev=mean(mean_sp_rev), mean_f1=mean(f1))

# All metric boxplots
homstrad.other_metrics <- homstrad.data %>%
  filter(type %in% c("irmsd", "nirmsd", "lddt", "apdb", "tc", "cs", "sp_fwd", "sp_rev")) %>%
  left_join(homstrad.counts, by=join_by(family)) %>% 
  filter(count >= 2) %>%
  mutate(score=as.numeric(score))
homstrad.other_metrics$tool <- factor(homstrad.other_metrics$tool, levels=unique(homstrad.other_metrics$tool[order(homstrad.other_metrics$base, homstrad.other_metrics$tool)]))

ggplot(
  homstrad.other_metrics
) +
  aes(y=score, x=tool, fill=base) +
  facet_wrap(~type, scales = "free") +
  geom_boxplot() +
  common_theme +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1))

ggsave(file=paste(BASEDIR, "figures/homstrad_all_metrics.pdf", sep=""), units="mm", width=300, height=300, dpi=300, bg="white")

# LDDT vs SoP/TC/CS reference-based scores on HOMSTRAD alignments
homstrad.wide <- homstrad.data %>%
  filter(tool != "HOMSTRAD" & type %in% c("lddt", "sp_fwd", "cs", "tc", "irmsd", "nirmsd", "apdb")) %>%
  pivot_wider(names_from=type, values_from=score, values_fn=as.numeric) %>%
  left_join(homstrad.counts, by=c("family")) %>%
  filter(count >= 4)

homstrad.longer <- homstrad.wide %>%
  pivot_longer(cols=c(sp_fwd, cs, tc, apdb, irmsd, nirmsd), names_to="score_type", values_to="score") %>%
  mutate(score_type = factor(score_type, levels=c('sp_fwd', 'tc', 'cs', 'apdb', 'irmsd', 'nirmsd')))

homstrad.correlations <- homstrad.longer %>% 
  group_by(score_type) %>%
  summarize(correlation=cor(x=score, y=lddt))
  
ggplot(homstrad.longer %>% filter(count >=4)) +
  aes(x=lddt, y=score, colour=score_type) + 
  facet_wrap(
    ~score_type,
    nrow=6,
    strip.position = "left", 
    labeller = as_labeller(c(
      sp_fwd = "Sum of pairs (SoP)",
      tc = "Total columns (TC)",
      cs = "Column score (CS)",
      apdb = "APDB",
      irmsd = "iRMSD (A)",
      nirmsd = "niRMSD (A)"
    )),
    scales="free_y"
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth=0.5, color="black") +
  geom_point(size=0.4, alpha=.5) +
  geom_text(
    data = homstrad.correlations,
    aes(x=Inf, y=Inf, label = paste("R =", round(correlation, 2))),
    color="black",
    hjust = 1.1,
    vjust = 5.0,
    size = 2,
  ) + 
  labs(
    title="Correlation of scoring metrics vs LDDT on Homstrad families with >=4 proteins",
    x="LDDT",
    y=element_blank()
  ) +
  guides(shape="none", colour="none") +
  common_theme +
  theme(
    strip.background=element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=6)
  )

ggsave(file=paste(BASEDIR, "figures/homstrad_correlations.pdf", sep=""), units="mm", width=180, height=170, dpi=300, bg="white")


# Panel 2: 1000 AFDB Clusters

# Read in long form score data
afdb = new.env()
afdb.data_base <- read.delim(
  paste(BASEDIR, "data/afdb_scores.tsv", sep=""),
  sep="\t",
  header=FALSE,
  col.names=c("family", "tool", "type", "score")
) %>%
  filter(type != "lddt_columns") %>%
  mutate(
    base=ifelse(name_map[tool] %in% structure_tools, "Structure-based", "Sequence-based"),
    score=as.numeric(score)
    # score=ifelse(type != "lddt_columns", as.numeric(score), score)
  )

# Cluster statistics
# need num_domains
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
  filter(type == "lddt") %>%
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
  filter(tool %in% c("foldmason", "foldmason_refine100") & type == "lddt") %>%
  dplyr::select(family, tool, score) %>%
  pivot_wider(names_from=tool, values_from=score) %>%
  mutate(delta = `foldmason_refine100` - `foldmason`)
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
afdb.summary <- summariseData(afdb.data_base, afdb.summary_time)
afdb.summary %>%
  dplyr::select(tool, avg_lddt, avg_lddt_1dom, avg_lddt_2dom, avg_lddt_3dom, avg_lddt_4dom, avg_time)


ggplot(afdb.data_base %>% filter(type == "lddt") %>% mutate(tool=name_map[tool])) +
  aes(x=reorder(tool, score, FUN=median), color=tool, y=score) +
  geom_boxplot() +
  labs(x="Tool", y="MSA LDDT") +
  scale_x_discrete(guide=guide_axis(angle=90)) +
  guides(color="none") +
  common_theme


x_max=15000
ymin=0
ymax=.8
afdb.plot_a <- ggplot(afdb.summary) +
  aes(x=speedup, y=avg_lddt, color=tool, shape=base) +
  geom_linerange(
    aes(ymin=first_quartile, ymax=third_quartile, color=NULL),
    linewidth=0.2,
    position=position_dodge(0.05),
    alpha=0.5
  ) +
  geom_point(size=1.2) +
  scale_x_log10(limits=c(0.5, x_max), expand=c(0,0), labels=label_log()) +
  scale_size_manual(values=c(2, 2.6)) +
  scale_y_continuous(limits=c(ymin, ymax), expand=c(0,0)) +
  scale_color_manual(values = colour_map) +
  scale_shape_manual(
    values = c("Sequence-based" = 16, "Structure-based" = 17),
    labels = c("Sequence-based" = "Sequence-based", "structure" = "Structure-based"),
  ) +
  labs(
    x="Speedup",
    y="LDDT",
    color="Tool",
    label=NULL,
    shape="Input",
    tag="ii)"
  ) +
  guides(
    shape=guide_legend(nrow=2),
    color=guide_legend(
      nrow=2,
      override.aes = list(shape = shape_map[afdb.summary$base], linetype = 0)
    ),
  ) +
  common_theme

afdb.plot_a

# Calculate per-tool scores wrt number of domains
afdb.plot_b <- ggplot(
  afdb.data_base %>%
    filter(type == 'lddt') %>%
    group_by(tool, num_domains) %>%
    summarise(avg_score = mean(score), base=first(base), .groups="keep") %>%
    mutate(tool=name_map[tool])
  ) +
  aes(x=num_domains, y=avg_score, color=tool, shape=base, group=tool) +
  scale_x_discrete(expand=c(.1,.1)) +
  scale_y_continuous(limits=c(ymin, ymax), expand=c(0,0)) +
  scale_color_manual(values = colour_map) +
  scale_shape_manual(
    values = c("Sequence-based" = 1, "Structure-based" = 2),
    labels = c("Sequence-based" = "Sequence-based", "Structure-based" = "Structure-based"),
  ) +
  geom_point(size=0.8, stroke=0.3) +
  geom_line(linetype='dashed', alpha=0.5, linewidth=0.2) +
  labs(
    x="Number of Domains",
    y="LDDT",
    color="Tool",
    shape="Input"
  ) +
  guides(shape="none", color="none") +
  common_theme

afdb.plot_a + afdb.plot_b + plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave(file=paste(BASEDIR, "newlddt.pdf", sep=""), units="mm", width=160, height=80, dpi=300)


# Breakdown of all tool LDDT scores on all AFDB clusters
family_order <- afdb.data_base %>%
  filter(type == "lddt", tool == "foldmason") %>%
  group_by(family) %>%
  summarise(ordering_score = mean(score), .groups = "drop") %>%
  arrange(desc(ordering_score)) %>%
  pull(family)

ggplot(afdb.data_base %>% filter(type == "lddt") %>% group_by(family, tool, num_domains, base) %>% summarise(avg=mean(score), .groups='drop') %>% group_by(num_domains) %>% mutate(tool=name_map[tool], family = factor(family, levels = family_order))) +
  facet_wrap(vars(num_domains), dir='v', ncol=1, scales='free_y') +
  scale_y_discrete(drop = TRUE) +
  aes(y=family, x=avg, color=tool, group=tool, shape=base) +
  #aes(y=reorder(family, avg, FUN=median), x=avg, color=tool, group=tool, shape=base) +
  scale_color_manual(name="tool", values = colour_map) +
  #geom_line(linetype=2, linewidth=0.1) +
  geom_point(size=3) +
  theme_minimal()

ggsave(file=paste(BASEDIR, "figures/afdb_all_families.pdf", sep=""), units="mm", limitsize=F, width=300, height=1000, bg="white")


# Per-family FoldMason refinement delta 
fmonly <- afdb.data_base %>% filter(type == 'lddt', tool=='foldmason')
fmronly <- afdb.data_base %>% filter(type == 'lddt', tool=='foldmason_refine100')
fmonly %>% inner_join(fmronly, by=c('family')) %>%
  mutate(diff=score.y-score.x) %>%
  dplyr::select(family, score.x, score.y, diff) %>%
  arrange(by=diff)


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

scientific_10 <- function(x) {
  xx <- dplyr::case_when(x>=0 & x<100 ~ as.character(x),
                         x>=100 & log10(x)%%1==0 ~ gsub(".+e\\+", "10^", scales::scientific_format()(x)),
                         TRUE ~ gsub("e\\+", " %*% 10^", scales::scientific_format()(x))) 
  parse(text = xx)
}

speed.times %>%
  dplyr::select(subset, tool, time) %>%
  pivot_wider(names_from=tool, values_from=time) %>%
  dplyr::select(subset, `3D-Coffee`, FAMSA, FoldMason)

options(scipen=999)
speed_plot <- ggplot(speed.times %>% filter(subsetNum %in% c(10, 100, 1000, 10000, 100000))) +
  aes(x=subsetNum, y=time, color=tool, group=tool, shape=base) +
  geom_point(size=1.2) +
  geom_line(alpha=0.5, linetype='dotted') +
  scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000), labels=scientific_10) +
  scale_y_log10(breaks=c(1, 10, 60, 600, 3600, 14400, 36000, 86400), labels=convert_time) +
  scale_color_manual(name="tool", values = colour_map) +
  scale_shape_manual(name="tool", values = shape_map) +
  labs(x="Subset size", y="Time (s)", color="Tool", tag="iii)") +
  guides(shape="none", colour="none") +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

speed_plot


# Plot the main figure from manuscript
design <- "
  1223
"
homstrad.plot + (afdb.plot_a + afdb.plot_b + plot_layout(axes="collect_y", tag_level="new")) + speed_plot +
  plot_layout(guides="collect", design=design) & #, heights=c(1, 1), widths=c(1, 1)) +
  theme(legend.position="bottom")

# device=svg is required for correct import into figma, converts text to paths though
# svglite preserves text but has weirdness with path stroke/fill when imported into figma
ggsave(file=paste(BASEDIR, "figures/benchmarks.pdf", sep=""), units="mm", width=160, height=60, dpi=300, bg="white")
ggsave(file=paste(BASEDIR, "figures/benchmarks.svg", sep=""), device=svg, units="mm", width=160, height=60, dpi=300, bg="white")


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

ggsave(file=paste(BASEDIR, "figures/homstrad_distribution.pdf", sep=""), units="mm", width=100, height=80, dpi=300, bg="white")


# Total proteins in each database
sum(pdb_protein_size$count)
sum(homstrad_protein_ranges$count)

# Flexibility scores
flex = new.env()
flex.data = read.csv(paste(BASEDIR, "flexibility/scores.tsv", sep=""), sep="\t") %>%
  pivot_longer(cols=c(pair, msa), names_to="type", values_to="score") %>%
  pivot_wider(names_from=method, values_from=score)

flex.data

# AFDB other metrics
afdb.all_copy <- afdb.data_base %>%
  filter(type %in% c("irmsd", "nirmsd", "apdb", "lddt")) %>%
  mutate(num_domains = "All")
afdb.other_metrics <- afdb.data_base %>%
  filter(type %in% c("irmsd", "nirmsd", "apdb", "lddt")) %>%
  bind_rows(afdb.all_copy) %>%
  mutate(num_domains = factor(num_domains, levels=c("All", 1, 2, 3, 4)))

afdb.other_metrics

# To display sequence tools then structure tools
afdb.other_metrics$tool <- factor(
  afdb.other_metrics$tool,
  levels=unique(afdb.other_metrics$tool[order(afdb.other_metrics$base, afdb.other_metrics$tool)])
)

ggplot(
  afdb.other_metrics %>% filter(type != "nirmsd" | score < 4.0)
) +
  aes(y=score, x=tool, fill=base) +
  facet_grid(type~num_domains, scales = "free_y", space = "fixed") +
  geom_boxplot() +
  common_theme +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1))

ggsave(file=paste(BASEDIR, "afdb_allmetrics_newlddt.pdf", sep=""), units="mm", width=180, height=160, dpi=300, bg="white")


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
  # ylim(0, 2500) +
  geom_boxplot() +
  facet_grid(rows = vars(num_domains), scales="free_y", space="fixed") +
  common_theme

ggsave(file=paste(BASEDIR, "figures/afdb_length_boxplots.pdf", sep=""), units="mm", width=180, height=160, dpi=300, bg="white")


# Identity distributions of AFDB cluster dataset
identities <- read.delim(
  paste(BASEDIR, "data/afdb_identities.csv", sep=""),
  sep=",",
  header=FALSE,
  col.names=c("family", "tool", "sequence_a", "sequence_b", "score")
)

ggplot(identities) +
  aes(x=score, y=after_stat(density), color=tool) +
  geom_histogram(fill='white', binwidth=0.01)  +
  labs(x="Sequence (MMseqs2) and structure (Foldseek) identity (%)", y="Frequency", color="Tool") +
  scale_color_discrete(labels=c('foldseek' = 'Foldseek', 'mmseqs' = 'MMseqs2')) +
  common_theme
  
ggsave(file=paste(BASEDIR, "figures/afdb_identities.pdf", sep=""), units="mm", width=180, height=160, dpi=300, bg="white")


# Homstrad pairwise F1 scores
pairs = new.env()
pairs.data <- read.csv(
  paste(BASEDIR, "data/homstrad_pair_scores.tsv", sep=""),
  sep='\t',
  col.names = c('family', 'tool', 'score'),
  header = F
)

pairs.fm_scores <- pairs.data %>%
  filter(tool == 'foldmason') %>%
  rename(fm_score = score) %>%
  dplyr::select(family, fm_score)

pairs.data <- pairs.data %>%
  filter(tool != "foldmason") %>%
  right_join(pairs.fm_scores, by="family")

pairs.counts <- pairs.data %>%
  group_by(tool) %>%
  summarise(N = sum(!is.na(score)), corr=cor(fm_score, score, use="complete.obs") %>% round(2))

ggplot(df2) +
  aes(x=fm_score, y=score, color=tool) +
  geom_point(size=0.8, alpha=0.5, stroke=NA) +
  geom_abline(color="black", alpha=0.5, linetype=2) +
  geom_text(data=counts, color="black", size=2, aes(x=.21, y=.96, label=paste0("N = ", N, "\nR = ", corr))) +
  facet_wrap(~ tool, labeller = as_labeller(function(tool) name_map[tool])) +
  labs(x="FoldMason F1 score", y="Tool F1 score") +
  guides(color="none") +
  common_theme

ggsave(file=paste(BASEDIR, "figures/homstrad_pair_f1_scores.pdf", sep=""), units = "mm", width = 180, height = 160)
