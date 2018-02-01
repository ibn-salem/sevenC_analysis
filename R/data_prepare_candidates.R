#*******************************************************************************
# Read CTCF motif sites from JASPAR TSV files and provide it as GRanges object
#*******************************************************************************

library(sevenC)  # to import() BED files
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
require(BSgenome.Hsapiens.UCSC.hg19) # for human genome sequence
require(DNAshapeR)    # for DNA shape prediction
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()

#*******************************************************************************
# Parameters ----
#*******************************************************************************
MOTIF_PVAL <- 2.5 * 1e-06
SHAPE_WINDOW = 500

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
  LoopRao2014_GM12878_File, seqinfo = seqInfoHg19)

trueLoopsTang2015 <- do.call(
  "c",
  map(LoopTang2015_GM12878_Files, 
         sevenC::parseLoopsTang2015, 
         seqinfo = seqInfoHg19))

gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
gi$loop <- factor(
  gi$Loop_Tang2015_GM12878 == "Loop" | gi$Loop_Rao_GM12878 == "Loop",
  c(FALSE, TRUE),
  c("No loop", "Loop")
)
# save file for faster reload
write_rds(gi, paste0(dataCandidatesPreifx, ".gi.rds"))

#*******************************************************************************
# Output motif sites as BED for loop and non-loop sites
#*******************************************************************************
ancGR <- regions(gi)
loops <- gi[gi$loop == "Loop"]
loopAncIdx <- unique(c(
  anchors(loops, "first", id = TRUE),
  anchors(loops, "second", id = TRUE)
))

loopAncGR <- ancGR[loopAncIdx]
nonLoopAncGR <- ancGR[-loopAncIdx]

export.bed(loopAncGR, paste0(dataCandidatesPreifx, ".loop_anchors.bed"))
export.bed(nonLoopAncGR, paste0(dataCandidatesPreifx, ".non-loop_anchors.bed"))

