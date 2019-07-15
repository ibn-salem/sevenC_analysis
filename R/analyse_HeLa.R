################################################################################
# Predict chromatin loopsin HeLa cells and validate with Hi-C and ChIA-PET
################################################################################


require(sevenC)    # devtools::install_github("ibn-salem/sevenC")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(rtracklayer)  # to import() BED files
require(rsample)
require(pryr) # for object_size()
require(feather)      # for efficient storing of data.frames
library(Vennerable)   # for Venn-Diagram

# require(multidplyr)   # for partition() and collect() to work in parallel

source("R/sevenC.functions.R")


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(10, parallel::detectCores() - 1)

# MIN_MOTIF_SIG <- 5
MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", paste0("v05_HeLa.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)



# v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1
screenPrefix <- file.path("results", paste0("v05_screen_TF_lfc.", 
                                            paste0("motifPval", MOTIF_PVAL), 
                                            "_w", WINDOW_SIZE, 
                                            "_b", BIN_SIZE))


# True loops in HeLa from Rao et al:
LoopRao2014_HeLa_File <- 
  "data/Rao2014/GSE63525_HeLa_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_HeLa_Files <- c(
  "data/Tang2015/GSM1872888_HeLa_CTCF_PET_clusters.txt",
  "data/Tang2015/GSM1872889_HeLa_RNAPII_PET_clusters.txt")



# metadata file
metaFile <- "data/ENCODE/metadata.fc_HELA_selected.tsv"

SELECTED_TF <- c(
  "RAD21",
  "CTCF",
  "ZNF143",
  "STAT1",
  "EP300",
  "POLR2A"
)

COL_SELECTED_TF = brewer.pal(length(SELECTED_TF), "Set1")
COL_SELECTED_TF_1 = brewer.pal(12, "Paired")[c(1, 3, 5, 7, 9, 11)]
COL_SELECTED_TF_2 = brewer.pal(12, "Paired")[c(2, 4, 6, 8, 10, 12)]
names(COL_SELECTED_TF_1) <- SELECTED_TF
names(COL_SELECTED_TF_2) <- SELECTED_TF

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

# parse ucscMeta file
meta <- read_tsv(metaFile,
                 col_types = 
                   cols(
                     `File accession` = col_character(),
                     TF = col_character(),
                     `Output type` = col_character(),
                     file_nrep = col_character(),
                     exp_nrep = col_integer(),
                     Lab = col_character(),
                     filePath = col_character()
                   )
)

# reformat metadata
meta <- meta %>% 
  filter(TF %in% SELECTED_TF) %>% 
  # mutate(name = paste0(TF, "_lfc")) %>%
  mutate(name = TF) %>%
  select(TF, name, filePath, everything())

# Output table with accession numbers for all TFs ------------------------------
access_table <- meta %>%
  select(TF, `Biosample term name`, `Output type`, `File accession`, `File download URL`, `Experiment accession`, everything()) %>% 
  select(-usedTF, -filePath, -file_nrep, -exp_nrep, -output_type, -name, -size) %>% 
  write_tsv(paste0(outPrefix, ".access_table.tsv"))

# Select motifs and parse input data -----------------------------------

# read CTCF moitf pairs as candidates
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))

# remove annotation from GM12878 cells
gi$loop_GM12878 <- gi$loop
mcols(gi)[, c("Loop_Rao_GM12878", "Loop_Tang2015_GM12878_CTCF", "Loop_Tang2015_GM12878_RNAPII", "loop")] <- NULL


# parse true loops in HeLa
seqInfo <- seqinfo(gi)
trueLoopsRao <- sevenC::parseLoopsRao(
  LoopRao2014_HeLa_File, seqinfo = seqInfo)
trueLoopsTang2015 <- do.call(
  "c",
  lapply(LoopTang2015_HeLa_Files, 
         sevenC::parseLoopsTang, 
         seqinfo = seqInfo))

gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_HeLa")
gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_HeLa")
gi$loop <- factor(
  gi$Loop_Tang2015_HeLa == "Loop" | gi$Loop_Rao_HeLa == "Loop",
  c(FALSE, TRUE),
  c("No loop", "Loop")
)

# save file for faster reload
write_rds(gi, paste0(outPrefix, ".gi_raw.rds"))
# gi <- read_rds(paste0(outPrefix, ".gi_raw.rds"))

# Count positive loops  --------------------------------------------------------
countLoopsDF <- gi %>% 
  mcols() %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  count(loop) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  write_tsv(paste0(outPrefix, ".count_positive_loops.tsv"))


# Get performance of GM12878 loops in HeLa -------------------------------------
conf_tab <- table(gi$loop, gi$loop_GM12878)
GM_sens <- sum(gi$loop == "Loop" & gi$loop_GM12878 == "Loop") / sum(gi$loop == "Loop")
GM_ppv <- sum(gi$loop == "Loop" & gi$loop_GM12878 == "Loop") / sum(gi$loop_GM12878 == "Loop")
GM_spec <- sum(gi$loop == "No loop" & gi$loop_GM12878 == "No loop") / sum(gi$loop == "No loop")

gm_performance <- tibble(
  name = "GM12878_measured",
  sens = GM_sens,
  ppv = GM_ppv,
  spec = GM_spec
)

write_tsv(gm_performance, str_c(outPrefix, ".gm_performance.tsv"))

# overlap of GM12878 loops and HeLa loops
loop_list <- list(
  GM12878 = which(gi$loop_GM12878 == "Loop"),
  HeLa = which(gi$loop == "Loop")
)
v <- Venn(loop_list)
pdf(str_c(outPrefix, ".loop_overlap_GM12878_HeLa.pdf"))
  plot(v, show = list(Faces = FALSE))
dev.off()

# Annotae with coverage and correlation ----------------------------------------

if (!GI_LOCAL ) {
  
  # iterate over all ChIP-seq sata sets
  for (i in seq_len(nrow(meta))) {
    
    message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
    
    #add coverage and correlation of coverage
    gi <- addCor(
      gi,
      meta$filePath[[i]],
      meta$name[[i]],
    )
  }  
  # Annotae with correlation across TFs -------------------------------------
  
  # get vector with coverage in whole anchor regions
  covCols <- meta$name
  covDF <- as_tibble(as.data.frame(mcols(regions(gi))[,covCols])) %>% 
    dplyr::mutate_all(.funs = function(l) map_dbl(l, sum))
  
  mcols(regions(gi))[, "cov_sum"] <- NumericList(as_tibble(t(covDF)))
  
  gi <- addCovCor(
    gi, 
    datacol = "cov_sum",
    colname = "across_TFs"
  )
  
  # save file for faster reload
  write_rds(gi, paste0(outPrefix, ".gi.rds"))

  df <- as_tibble(as.data.frame(mcols(gi))) %>%
    mutate(
      id = 1:nrow(.)
    ) %>% 
    select(id, loop, everything()) 
  
  write_feather(df, paste0(outPrefix, ".df.feather"))
  
} else {
  gi <- read_rds(paste0(outPrefix, ".gi.rds"))
  df <- read_feather(paste0(outPrefix, ".df.feather"))
}

#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------

#===============================================================================
# Number of positive looops in HeLA data
#===============================================================================
countDF <- df %>% 
  count(loop) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  write_tsv(paste0(outPrefix, ".number_of_loops.tsv"))

#===============================================================================
# Training and prediction in cross-validation
#===============================================================================

# get design formula for each TF and model
designDF <- tibble(
  name = factor(meta$name, meta$name),
  design = map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x))),
  color = COL_SELECTED_TF_2
  )

write_rds(designDF, paste0(outPrefix, "designDF.rds"))
# designDF <- read_rds(paste0(outPrefix, "designDF.rds"))

TFspecific_ModelDF <- read_tsv(paste0(screenPrefix, ".TFspecific_ModelDF.tsv"))
allTfModelDF <- read_tsv(paste0(screenPrefix, ".allTfModelDF.tsv"))
bestNModelDF <- read_tsv(paste0(screenPrefix, ".bestNModelDF.tsv"))

# pie(rep(1, nrow(designDF)), col = designDF$color, labels = designDF$name)

predDF <- designDF %>% 
  # add TF specific model as colum list
  mutate(
    betas_specificTF = map(name, ~ pull(filter(TFspecific_ModelDF, TF == .x), estimate_mean)),
  ) %>% 
  # add predictions using specific TF and other models as columns
  mutate(
    pred_specificTF = map2(design, betas_specificTF, 
                           ~ sevenC:::predLogit(df, formula = .x, betas = .y)),
    pred_allTF = map(design, ~ sevenC:::predLogit(df, .x, allTfModelDF$estimate_mean)),
    pred_bestN = map(design, ~ sevenC:::predLogit(df, .x, bestNModelDF$estimate_mean))
  )
  
write_rds(predDF, path = paste0(outPrefix, "predDF.rds"))
# predDF <- read_rds(paste0(outPrefix, "predDF.rds"))

# combine pred types
evalDF <- predDF %>% 
  gather(key = "pred_type", value = "pred", starts_with("pred_")) %>% 
  mutate(
    modnames = paste(name, str_replace(pred_type, "pred_", ""), sep = "_")
  )

#===============================================================================
# Performance Evaluation
#===============================================================================

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = evalDF$pred,
  labels = list(df$loop),
  modnames = evalDF$modnames,
  posclass = levels(df$loop[[1]])[2],
  x_bins = 100)

# prc_base = precrec:::.get_pn_info(curves)$prc_base
prc_base = mean(df$loop == "Loop")
write_rds(prc_base, paste0(outPrefix, "prc_base.rds"))
# prc_base <- read_rds(paste0(outPrefix, "prc_base.rds"))

# get data.frame with auc values
aucDF <-  as_tibble(auc(curves)) %>% 
  separate(modnames, c("name", "pred_type"), sep = "_") %>% 
  mutate(name = factor(name, SELECTED_TF)) %>% 
  arrange(pred_type, name)


write_feather(aucDF, paste0(outPrefix, ".aucDF.feather"))
# aucDF <- read_feather(paste0(outPrefix, ".aucDF.feather"))

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
plotDF <- aucDF %>% 
  mutate(pred_type = factor(
    pred_type, 
    levels = c("specificTF", "allTF", "bestN"),
    labels = c("Specific TF", "Avg. all TF", "Avg. best 10 TF")),
  curvetypes = factor(curvetypes, c("ROC", "PRC"))
  )

p <- ggplot(plotDF, aes(x = name, y = aucs, fill = name)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  # geom_text(aes(label = round(aucs, 2), y = aucs), size = 3, vjust = 1.2, angle = 0) +
  geom_text(aes(label = round(aucs, 2), y = aucs), angle = 90, hjust = 1.2) +
  facet_grid(curvetypes ~ pred_type, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_SELECTED_TF_2) +
  labs(x = "ChIP-seq data", y = "Prediction performance (AUC)")

ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC_by_predtype.barplot.pdf"), w = 5, h = 5)


aucDFspecific <- aucDF %>%
  filter(pred_type == "specificTF", curvetypes == "PRC")

p <- ggplot(aucDFspecific, aes(x = name, y = aucs, fill = name)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = round(aucs, 2), y = aucs), size = 5, vjust = 1.2, angle = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = COL_SELECTED_TF_2) +
  labs(x = "", y = "Prediction performance\n(AUC PRC)")
ggsave(p, file = paste0(outPrefix, ".AUC_PRC_specificTF.barplot.pdf"), w = 3.5, h = 7)

#-------------------------------------------------------------------------------
# get ROC and PRC plots
#-------------------------------------------------------------------------------
curveDF <- as_tibble(as.data.frame(curves)) %>% 
  separate(modname, c("name", "pred_type"), sep = "_") %>% 
  mutate(name = factor(name, SELECTED_TF)) 

write_feather(curveDF, paste0(outPrefix, ".curveDF.feather"))
# curveDF <- read_feather(paste0(outPrefix, ".curveDF.feather"))

# subsample curve data
curveDFsub <- curveDF %>%
  group_by(pred_type, name, type) %>% 
  sample_n(1000) %>% 
  ungroup()
  
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC", pred_type == "specificTF")

rocDF <- curveDFsub %>% 
  filter(type == "ROC", pred_type == "specificTF")

g <- ggplot(rocDF, aes(x = x, y = y, color = name)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dotted") +
  geom_point(aes(x = 1 - spec, y = sens), data = gm_performance, size = 2, color = "darkgray") +
  geom_text(aes(x = 1 - spec, y = sens, label = "GM12878 loops"), data = gm_performance, 
            color = "darkgray", hjust = 1.05, vjust = -0.5, angle = 90) +
  theme_bw() +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual("",
                     values = COL_SELECTED_TF_2,
                     labels = paste0(
                       aucDFroc$name, 
                       ": AUC=", 
                       signif(aucDFroc$aucs, 3)
                     ),
                     guide = guide_legend(
                       override.aes = list(size = 2))
  ) +
  theme(aspect.ratio = 1, text = element_text(size = 15), legend.position = c(.7,.3),
        legend.background = element_rect(fill = alpha('white', 0.1)))
ggsave(g, file = paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

#-------------------------------------------------------------------------------
# get PRC plots
#-------------------------------------------------------------------------------
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC", pred_type == "specificTF")

prc_base = precrec:::.get_pn_info(curves)$prc_base
write_rds(prc_base, paste0(outPrefix, "prc_base.rds"))
# prc_base <- read_rds(paste0(outPrefix, "prc_base.rds"))

prcDF <- curveDFsub %>% 
  filter(type == "PRC", pred_type == "specificTF")

g <- ggplot(prcDF, aes(x = x, y = y, color = name)) +
  geom_line() +
  geom_abline(intercept = prc_base, slope = 0, lty = "dotted") +
  geom_point(aes(x = sens, y = ppv), data = gm_performance, size = 2, color = "darkgray") +
  geom_text(aes(x = sens, y = ppv, label = "GM12878 loops"), data = gm_performance, 
            color = "darkgray", hjust = 1.1) +
  theme_bw() +
  labs(x = "Recall", y = "Precision") +
  scale_color_manual("",
                     values = COL_SELECTED_TF_2,
                     labels = paste0(
                       aucDFroc$name, 
                       ": AUC=", 
                       signif(aucDFroc$aucs, 3)
                     ),
                     guide = guide_legend(
                       override.aes = list(size = 2))
  ) +
  theme(aspect.ratio = 1, 
        legend.position = "none", 
        text = element_text(size = 15)) +
  ylim(0, 1)
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)


#*******************************************************************************
# Binary prediction ------
#*******************************************************************************

# bestNModelDF <- read_tsv(paste0(screenPrefix, ".bestNModelDF.tsv"))

# read TF specific models
TFspecific_ModelDF <- read_tsv(paste0(screenPrefix, ".TFspecific_ModelDF.tsv"))

Rad21_model <- TFspecific_ModelDF %>% 
  filter(TF == "RAD21")

BestN_cutof <- read_tsv(paste0(screenPrefix, ".topNf1ModelDF.tsv")) %>% 
  pull(mean_max_cutoff)


# add predictions using RAD21 model
gi <- predLoops(
  gi,
  formula = loop ~ dist + strandOrientation + score_min + cor_RAD21,
  betas = Rad21_model$estimate_mean,
  colname = "pred_Rad21",
  cutoff = NULL
)

# add binary predictions
mcols(gi)$predBinary_Rad21 <- !is.na(mcols(gi)$pred_Rad21) & mcols(gi)$pred_Rad21 >= BestN_cutof


#*******************************************************************************
# Write predicted interactions as genome browser track
#*******************************************************************************


# write only subset of interacting pairs
loopGI <- gi[gi$loop == "Loop"]
writeLongRangeTrackFormat(
  gi = loopGI, 
  score_vec = rep(1, length(loopGI)), 
  output_file = paste0(outPrefix, ".gi.HeLa.loop_sub.longrange_track.bed")
)

# write only subset of predicted loops
rad21GI <- gi[gi$predBinary_Rad21]
writeLongRangeTrackFormat(
  gi = rad21GI, 
  score_vec = 100 * rad21GI$pred_Rad21, 
  output_file = paste0(outPrefix, ".gi.HeLa.pred_Rad21_sub.longrange_track.bed")
)

#===============================================================================
#===============================================================================
