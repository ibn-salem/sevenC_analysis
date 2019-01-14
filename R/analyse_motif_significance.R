#*******************************************************************************
# Anaylse the motif significance threshold of CTCF moitfs from JASPAR
#*******************************************************************************


library(sevenC)  # to import() BED files
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()
require(RColorBrewer)   # for nice colors

#*******************************************************************************
# Parameters ----
#*******************************************************************************
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

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
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", "CTCF_JASPAR.v01.motif_sig")

#*******************************************************************************
# Parse CTCF motif sites from JASPAR track -----------------
#*******************************************************************************

# define colnames and parse TSV file from JASPAR
# header: chr `start (1-based)`   end `rel_score * 1000` `-1 * log10(p_value) * 100` strand
col_names = c("chr", "start", "end", "name", "score", "log10_pval_times_100", "strand")
allMotifDF <- read_tsv(JASPAR_HG19_CTCF, col_names = col_names, skip = 1, 
                    col_type = cols(
                      chr = col_character(),
                      start = col_integer(),
                      end = col_integer(),
                      name = col_character(),
                      score = col_integer(),
                      log10_pval_times_100 = col_integer(),
                      strand = col_character()
                    )) %>% 
  mutate(log10_pval = log10_pval_times_100 / 100)


#*******************************************************************************
# Analze the cutof threhold on p-value
#*******************************************************************************
df <- allMotifDF %>% 
  filter(log10_pval >= 5)

# get seqinfo object
seqInfoHg19 <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

# build GRanges object
motifGR <- GRanges(df$chr, IRanges(df$start, df$end),
                   strand = df$strand,
                   score = df$log10_pval,
                   seqinfo = seqInfoHg19)
# sort
motifGR <- sort(motifGR)


# parse loops
trueLoopsRao <- parseLoopsRao(LoopRao2014_GM12878_File, seqinfo = seqInfoHg19)

trueLoopsTang2015_list <- map(LoopTang2015_GM12878_Files, 
                              sevenC::parseLoopsTang, 
                              seqinfo = seqInfoHg19)

trueLoopsTang2015 <- do.call("c", trueLoopsTang2015_list)


gi <- prepareCisPairs(motifGR, maxDist = 1e+6)
gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
gi <- addInteractionSupport(gi, trueLoopsTang2015_list[[1]], "Loop_Tang2015_GM12878_CTCF")
gi <- addInteractionSupport(gi, trueLoopsTang2015_list[[2]], "Loop_Tang2015_GM12878_RNAPII")
# gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
gi$loop <- factor(
  gi$Loop_Rao_GM12878 == "Loop" | 
    gi$Loop_Tang2015_GM12878_CTCF == "Loop" | 
    gi$Loop_Tang2015_GM12878_RNAPII == "Loop",
  c(FALSE, TRUE),
  c("No loop", "Loop")
)


moitf_log10_pval <- motifGR$score
gi_max_score = gi$score_min
loop = gi$loop == "Loop"

motif_cutoff_DF <- tibble(
  p = 1e-6 * seq(1, 10, 0.5),
  log10_p = -log10(p),
  n_motif = map_int(log10_p, ~ sum(moitf_log10_pval >= .x)),
  n_gi =  map_int(log10_p, ~ sum(gi$score_min >= .x)),
  n_gi_loop = map_int(log10_p, ~ sum(loop[gi$score_min >= .x])),
  percent_loop = n_gi_loop / n_gi * 100
)

write_tsv(motif_cutoff_DF, paste0(outPrefix, ".moitf_cutoff_DF.tsv"))
# motif_cutoff_DF <- read_tsv(paste0(outPrefix, ".moitf_cutoff_DF.tsv"))

p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = n_motif)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -log10(MOTIF_PVAL), linetype = 2) +
  theme_bw() + 
  labs(x = "Motif significance\n-log10(p-value)", y = "Motifs")
ggsave(paste0(outPrefix, ".motifs_by_moitf_cutoff.pdf"), p,
       w = 2, h = 2)

p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = n_gi)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -log10(MOTIF_PVAL), linetype = 2) +
  theme_bw() + 
  labs(x = "Motif significance\n-log10(p-value)", y = "Motif pairs")
ggsave(paste0(outPrefix, ".motifs_pairs_by_moitf_cutoff.pdf"),
       w = 2, h = 2)

p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = percent_loop)) +
  geom_line(color = COL_LOOP[[2]]) +
  geom_point(color = COL_LOOP[[2]]) +
  geom_vline(xintercept = -log10(MOTIF_PVAL), linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = COL_LOOP[[2]]) +
  labs(x = "Motif significance\n-log10(p-value)", y = "Percent true loops")
ggsave(paste0(outPrefix, ".percent_true_loops_by_moitf_cutoff.pdf"),
       w = 2, h = 2)

