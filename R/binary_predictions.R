################################################################################
# Analysis of binary predictd chromatin looping interactions using 7C
################################################################################


library(sevenC)    # devtools::install_github("ibn-salem/sevenC")
library(rtracklayer)  # to import() BED files
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(modelr)       # for tidy modeling
library(precrec)      # for ROC and PRC curves
library(RColorBrewer)   # for nice colors
library(rsample)
library(pryr) # for object_size()
library(feather)      # for efficient storing of data.frames
library(ROCR)         # for binary clasification metrices

source("R/sevenC.functions.R")

# 0) Set parameter --------------------------------------------------------

MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

outPrefix <- file.path("results", paste0("v05_selected_models.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

screenPrefix <- file.path("results", paste0("v05_screen_TF_lfc.", 
                                            paste0("motifPval", MOTIF_PVAL), 
                                            "_w", WINDOW_SIZE, 
                                            "_b", BIN_SIZE))

# SELECTED_TF <- c(
#   "RAD21",
#   "CTCF",
#   "ZNF143",
#   "STAT1",
#   "EP300",
#   "POLR2A"
# )
# 
# COL_SELECTED_TF = brewer.pal(length(SELECTED_TF), "Set1")
# COL_SELECTED_TF_1 = brewer.pal(12, "Paired")[c(1, 3, 5, 7, 9, 11)]
# COL_SELECTED_TF_2 = brewer.pal(12, "Paired")[c(2, 4, 6, 8, 10, 12)]
# names(COL_SELECTED_TF_1) <- SELECTED_TF
# names(COL_SELECTED_TF_2) <- SELECTED_TF
# 
# 


#*******************************************************************************
# Binary prediction ------
#*******************************************************************************
gi <- read_rds(paste0(outPrefix, ".gi.rds"))

# read TF specific models
TFspecific_ModelDF <- read_tsv(paste0(screenPrefix, ".TFspecific_ModelDF.tsv"))

Rad21_model <- TFspecific_ModelDF %>% 
  filter(TF == "RAD21")

Rad21_cutof <- read_rds(paste0(screenPrefix, ".f1ModelDF.rds")) %>% 
  filter(TF == "RAD21") %>% 
  pull(max_cutoff)

# add predictions using RAD21 model
gi <- predLoops(
  gi,
  formula = loop ~ dist + strandOrientation + score_min + cor_RAD21,
  betas = Rad21_model$estimate_mean,
  colname = "pred_Rad21",
  cutoff = NULL
)

# add binary predictions
mcols(gi)$predBinary_Rad21 <- !is.na(mcols(gi)$pred_Rad21) & mcols(gi)$pred_Rad21 >= Rad21_cutof

#*******************************************************************************
# Write genome browser tracks   -----
#*******************************************************************************

# write CTCF moitfs as BED file
export.bed(regions(gi), paste0(outPrefix, ".motifs.bed"), index = TRUE)

onChr22 <- seqnames(regions(gi))[anchors(gi, type = "first", id = TRUE)] == "chr22"
chr22GI <- gi[onChr22]

# write all chr22 pairs (with labeld true ones)
writeLongRangeFormat(
  gi = chr22GI, 
  score_vec = ifelse(mcols(chr22GI)$loop == "Loop", 1, -1), 
  output_file = paste0(outPrefix, ".gi.loop.chr22.longrange.txt")
)

writeLongRangeTrackFormat(
  chr22GI,
  score_vec = ifelse(mcols(chr22GI)$loop == "Loop", 1, -1), 
  output_file = paste0(outPrefix, ".gi.loop.chr22.longrange_track.bed"),
  index = TRUE
)

# write all pairs to track format
writeLongRangeTrackFormat(
  gi,
  score_vec = ifelse(mcols(gi)$loop == "Loop", 1, -1), 
  output_file = paste0(outPrefix, ".gi.loop.longrange_track.bed")
)

# write all chr22 pairs (with labeld predictions)
writeLongRangeFormat(
  gi = chr22GI, 
  score_vec = ifelse(chr22GI$predBinary_Rad21, 1, -1), 
  output_file = paste0(outPrefix, ".gi.pred_Rad21.chr22.longrange.txt")
)

# write only subset of interacting pairs
loopGI <- gi[gi$loop == "Loop"]
writeLongRangeTrackFormat(
  gi = loopGI, 
  score_vec = rep(1, length(loopGI)), 
  output_file = paste0(outPrefix, ".gi.loop_sub.longrange_track.bed")
)

# write only subset of predicted loops
rad21GI <- gi[gi$predBinary_Rad21]
writeLongRangeTrackFormat(
  gi = rad21GI, 
  score_vec = 100 * rad21GI$pred_Rad21, 
  output_file = paste0(outPrefix, ".gi.pred_Rad21_sub.longrange_track.bed")
)

# write to JuiceBox 2D annotation format
writeJuiceboxFormat(rad21GI, resolution = 1000, 
                    output_file = paste0(outPrefix, ".gi.pred_Rad21_sub.2D_annot3.tsv"))

#===============================================================================
#===============================================================================
