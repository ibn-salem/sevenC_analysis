
require(GenomicRanges)
require(rtracklayer)
require(tidyverse)

MOTIF_PVAL <- 2.5 * 1e-06

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
          paste0("CTCF_JASPAR.v02.pval_", MOTIF_PVAL))

peak_file <- "data/ENCODE/Peaks/ENCFF710VEH.bed"
out_file <- "data/ENCODE/Peaks/ENCFF710VEH.bed.CTCF_motif.bed"

motifGR <- read_rds(paste0(dataCandidatesPreifx, "motifGR.rds"))

peak_df <- read_tsv(peak_file, col_names = FALSE)
peak_gr <- GRanges(peak_df$X1, 
                   IRanges(peak_df$X2, peak_df$X3),
                   "*",
                   score = peak_df$X7)

# compute ovlerlap of motifs with peaks
ol <- findOverlaps(motifGR, peak_gr)

# add peak score
motifGR$peak_score <- NA
motifGR$peak_score[queryHits(ol)] <- peak_gr$score[subjectHits(ol)]
# normalize to be in [0, 1000]
motifGR$peak_score <- motifGR$peak_score / max(motifGR$peak_score, na.rm = TRUE) * 1000

# use peak score as score
motifGR$score <- motifGR$peak_score

# filter for only motifs in peaks
motifGR <- motifGR[!is.na(motifGR$score)]

# write motifs with peak socre as BED file
export.bed(motifGR, out_file)
