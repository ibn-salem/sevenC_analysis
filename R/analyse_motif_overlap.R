#*******************************************************************************
# Anaylse overlap of measured loops with CTCF motif pairs
#*******************************************************************************

library(sevenC)  # to import() BED files
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()
require(RColorBrewer)   # for nice colors
library(UpSetR)

#*******************************************************************************
# Parameters ----
#*******************************************************************************
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

COL_LOOP_DATA = brewer.pal(8, "Set1")[c(4,5,7)]

MOTIF_PVAL <- 2.5 * 1e-06

JASPAR_HG19_CTCF <- "data/JASPAR2018/MA0139.1.tsv"

# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_GM12878_Files <- c(
  "data/Tang2015/GSM1872886_GM12878_CTCF_PET_clusters.txt",
  "data/Tang2015/GSM1872887_GM12878_RNAPII_PET_clusters.txt")

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v02.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", "CTCF_JASPAR.v02.motif_sig")

#*******************************************************************************
# Parse CTCF motif sites from JASPAR track -----------------
#*******************************************************************************

gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))

# parse loops ------------------------------------------------------------------

# get seqinfo object
seqInfoHg19 <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)


trueLoopsRao <- parseLoopsRao(LoopRao2014_GM12878_File, seqinfo = seqInfoHg19)

trueLoopsTang2015_list <- map(LoopTang2015_GM12878_Files, 
                              sevenC::parseLoopsTang, 
                              seqinfo = seqInfoHg19)

# Analyse overalp of true loops with motif pairs -------------------------------
dataset_names <- c("Rao2014 Hi-C", "Tang2015 CTCF ChIA-PET", "Tang2015 PolII ChIA-PET")


loop_data_sets <- tibble(
  name = factor(dataset_names, dataset_names),
  loops = list(trueLoopsRao, trueLoopsTang2015_list[[1]], trueLoopsTang2015_list[[2]]),
  loop_n = map_int(loops, length),
  overlaps_motif_pairs = map(loops, overlapsAny, gi),
  overlap_counts = map(loops, countOverlaps, gi),
  overlap_n = map_int(overlaps_motif_pairs, sum),
  overlap_percent = overlap_n / loop_n * 100
  )

ggplot(loop_data_sets, aes(x = fct_rev(name), y = overlap_percent, fill = name, 
                                label = str_c(overlap_n, "\n(", round(overlap_percent, 2), "%)")
                                )) +
  geom_bar(stat = "identity") +
  geom_text(hjust = "inward") +
  coord_flip() +
  scale_fill_manual(values = COL_LOOP_DATA) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Percent overlap with CTCF pairs",
       x = "")
ggsave(paste0(outPrefix, ".True_loops_overlap_with_CTCF_pairs.pdf"), w = 4.5, h = 3)

loops <- loop_data_sets %>% 
  mutate(
    dist = map(loops, pairdist, type = "mid")
  ) %>% 
  select(name, dist, overlaps_motif_pairs, overlap_counts) %>% 
  unnest(dist, overlaps_motif_pairs, overlap_counts) %>% 
  mutate(
    dist_range = cut(dist, 10 ** (1:9))
  )

# distance distribution of mesarued loops
ggplot(loops, aes(x = dist / 1000, stat(count), fill = fct_rev(name), color = fct_rev(name))) +
  geom_density(alpha = 0.25) + 
  geom_vline(xintercept = 10^3, color = "red", linetype = 2) +
  scale_x_log10() +
  scale_fill_manual("", values = rev(COL_LOOP_DATA)) +
  scale_color_manual("", values = rev(COL_LOOP_DATA)) +
  theme_bw() +
  theme(
    # legend.position = "bottom",
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 8)
  ) + 
  guides(fill = guide_legend(nrow = 3)) +
  labs(x = "Distance [kb]", y = "Interactions")
ggsave(paste0(outPrefix, ".True_loops_by_distance.pdf"), w = 4.5, h = 3)

# number of CTCF moitf pairs per loop
overlap_counts <- loops %>% 
  mutate(
    overlap_counts_group = ifelse(overlap_counts <= 5, as.character(overlap_counts), ">5") %>% 
      factor(c(as.character(1:5), ">5"))
  ) %>% 
  group_by(name) %>% 
  count(name, overlap_counts_group) %>% 
  filter(overlap_counts_group != "0") %>% 
  mutate(
    percent = n / sum(n) * 100
    )

ggplot(overlap_counts, aes(x = fct_rev(name), y = percent , 
                           label = str_c(n, "\n(", round(percent, 2), "%)"), 
                           fill = overlap_counts_group)) +
  geom_bar(stat = "identity") + 
  geom_text(position = position_stack(vjust = 0.5), vjust = "center", alpha = 1, 
            data = filter(overlap_counts, overlap_counts_group %in% c("1", "2")), size = 3) +
  # scale_fill_manual("", values = COL_LOOP_DATA) +
  # scale_color_manual("", values = COL_LOOP_DATA) +
  scale_fill_brewer("Number of \n overlapping\nmotif pairs", palette = "BuGn") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)
  ) + 
  labs(x = "", y = "Interactions [%]")
ggsave(paste0(outPrefix, ".True_loops_overlap_with_CTCF_pairs_counts.pdf"), w = 4.5, h = 3)

# Overlap ob truth data sets----------------------------------------------------

loop_set_df <- mcols(gi) %>% 
  as.data.frame() %>% 
  select(Loop_Rao_GM12878, Loop_Tang2015_GM12878_CTCF, Loop_Tang2015_GM12878_RNAPII) %>% 
  set_names(dataset_names) %>% 
  mutate_all(function(x) ifelse(x == "Loop", 1, 0))

pdf(paste0(outPrefix, ".True_loops_dataset_overlap.pdf"), w = 4.5, h = 3, onefile=FALSE)
  upset(loop_set_df, sets = rev(dataset_names), keep.order = TRUE)
dev.off()
