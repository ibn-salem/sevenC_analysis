#*******************************************************************************
# Read CTCF motif sites from JASPAR TSV files and provide it as GRanges object
#*******************************************************************************


library(chromloop)  # to import() BED files
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()

#*******************************************************************************
# Parameters ----
#*******************************************************************************
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

dir.create(dirname(dataCandidatesPreifx), showWarnings = FALSE)
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

motifDF <- allMotifDF %>%
  filter(log10_pval >= -log10(MOTIF_PVAL))

#*******************************************************************************
# Build GRanges object using hg19 seqinfo object  -----------------
#*******************************************************************************
seqInfoHg19 <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

# build GRanges object
motifGR <- GRanges(motifDF$chr, IRanges(motifDF$start, motifDF$end),
                   strand = motifDF$strand,
                   score = motifDF$log10_pval,
                   seqinfo = seqInfoHg19)
# sort
motifGR <- sort(motifGR)

#*******************************************************************************
# Save moitf GRanges as .rds file
#*******************************************************************************
write_rds(motifGR, paste0(dataCandidatesPreifx, "motifGR.rds"))

#*******************************************************************************
# Prepare motif pairs as candidates ----
#*******************************************************************************

# get all pairs within 1Mb
gi <- prepareCisPairs(motifGR, maxDist = 10^6)

#*******************************************************************************
# Add true loop annotations ----
#*******************************************************************************

# parse loops
trueLoopsRao <- parseLoopsRao(
  LoopRao2014_GM12878_File, seqinfo = seqInfo)

trueLoopsTang2015 <- do.call(
  "c",
  map(LoopTang2015_GM12878_Files, 
         chromloop::parseLoopsTang2015, 
         seqinfo = seqInfo))

gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
gi$loop <- factor(
  gi$Loop_Tang2015_GM12878 == "Loop" | gi$Loop_Rao_GM12878 == "Loop",
  c(FALSE, TRUE),
  c("No loop", "Loop")
)
# save file for faster reload
save(gi, file = paste0(dataCandidatesPreifx, ".gi.Rdata"))


#*******************************************************************************
# Analze the cutof threhold on p-value
#*******************************************************************************

df <- allMotifDF %>% 
  filter(log10_pval >= 5)

# build GRanges object
motifGR <- GRanges(df$chr, IRanges(df$start, df$end),
                   strand = df$strand,
                   score = df$log10_pval,
                   seqinfo = seqInfoHg19)
# sort
motifGR <- sort(motifGR)
gi <- prepareCisPairs(motifGR, maxDist = 1e+6)
gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
gi$loop <- factor(
  gi$Loop_Tang2015_GM12878 == "Loop" | gi$Loop_Rao_GM12878 == "Loop",
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

write_tsv(motif_cutoff_DF, paste0(dataCandidatesPreifx, ".moitf_cutoff_DF.tsv"))

p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = n_motif)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(x = "Motif significance -log10(p-value)", y = "Number of motifs")
ggsave(paste0(dataCandidatesPreifx, ".motifs_by_moitf_cutoff.pdf"),
       w = 3, h = 3)

p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = n_gi)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(x = "Motif significance -log10(p-value)", y = "Number of motif pairs")
ggsave(paste0(dataCandidatesPreifx, ".motifs_pairs_by_moitf_cutoff.pdf"),
       w = 3, h = 3)


p <- ggplot(motif_cutoff_DF, aes(x = log10_p, y = percent_loop)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  labs(x = "Motif significance -log10(p-value)", y = "Percent true loops")
ggsave(paste0(dataCandidatesPreifx, ".percent_true_loops_by_moitf_cutoff.pdf"),
       w = 3, h = 3)
