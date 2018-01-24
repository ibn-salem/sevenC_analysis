################################################################################
# Analysis of predictd chromatin looping interactions using the sevenC tool
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
library(multidplyr)   # for partition() and collect() to work in parallel
library(ROCR)         # for binary clasification metrices

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


outPrefix <- file.path("results", paste0("v05_selected_models.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

# v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1
screenPrefix <- file.path("results", paste0("v05_screen_TF_lfc.", 
                                            paste0("motifPval", MOTIF_PVAL), 
                                            "_w", WINDOW_SIZE, 
                                            "_b", BIN_SIZE))


# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_GM12878_Files <- c(
  "data/Tang2015/GSM1872886_GM12878_CTCF_PET_clusters.txt",
  "data/Tang2015/GSM1872887_GM12878_RNAPII_PET_clusters.txt")


# metadata file
metaFile <- "data/ENCODE/metadata.fcDF.tsv"

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
# setup cluster
#-------------------------------------------------------------------------------

# partion data for parallel processing
cluster <- create_cluster(N_CORES) %>%
  cluster_library(packages = c("sevenC", "tidyverse"))

# evaluate help function code on each cluster
cluster_eval(cluster, source("R/sevenC.functions.R"))


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
  arrange(factor(name, SELECTED_TF)) %>% 
  select(TF, name, filePath, everything())


write_tsv(meta, paste0(outPrefix, ".meta.tsv"))
# meta <- read_tsv(paste0(outPrefix, ".meta.tsv"))

# Select motifs and parse input data -----------------------------------

# read CTCF moitf pairs as candidates
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))

# Annotae with coverage and correlation -------------------------------------

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
    
    # source("R/old_code.R")
    # regions(gi) <- addCovToGR_OLD(regions(gi), meta$filePath[[i]], colname = meta$name[[i]])
    # regions(gi) <- addCovToGR(regions(gi), meta$filePath[[i]], colname = meta$name[[i]])
    # regions(gi) <- addCovToGR_asGR(regions(gi), meta$filePath[[i]], colname = meta$name[[i]])
    # gi <- addCovCor(gi, datacol = meta$name[[i]], colname = paste0("cor_", meta$name[[i]]))
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
  
} else {
  gi <- read_rds(paste0(outPrefix, ".gi.rds"))
}

#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------
df <- as_tibble(as.data.frame(mcols(gi))) %>%
  mutate(
    id = 1:nrow(.)
  ) %>% 
  select(id, loop, everything()) 

write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

#===============================================================================
# Training and prediction in cross-validation
#===============================================================================

# k fold cross validation with 1 repeats
set.seed(3579)

tidyCV <- df %>% 
  vfold_cv(V = K, repeats = 1) %>% 
  tidy()

write_feather(tidyCV, paste0(outPrefix, ".tidyCV.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))


designList <- list(
  "Dist" =  loop ~ dist,
  "Orientation" = loop ~ strandOrientation,
  "Motif" = loop ~ score_min,
  "Dist+Orientation+Motif" = loop ~ dist + strandOrientation + score_min,
  "across_TFs" = as.formula("loop ~ dist + strandOrientation + score_min + across_TFs"),
  "all_TF" = as.formula(paste0("loop ~ dist + strandOrientation + score_min + ", paste(paste0("cor_", meta$name), collapse = " + ")))
)

# get design formula for each TF and model
designDF <- tibble(
  name = c(
    names(designList),
    paste0(meta$name, "_only"),
    meta$name
  ),
  design = c(
    designList,
    map(meta$name, ~as.formula(paste0("loop ~ cor_", .x)) ),
    map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x)) )
  ),
  color = c(
    c("greenyellow", "gold2", "khaki", "brown4", "darkgray", "grey30"),
    COL_SELECTED_TF_1,
    COL_SELECTED_TF_2
  )
) %>% 
  mutate(name = factor(name, c(names(designList)[1:4], paste0(meta$name, "_only"), meta$name, "all_TF", "across_TFs"))) %>% 
  arrange(name)

write_rds(designDF, paste0(outPrefix, "designDF.rds"))

# pie(rep(1, nrow(designDF)), col = designDF$color, labels = designDF$name)

# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = designDF$name, Fold) %>% 
  # add design formular for each TF
  left_join(designDF, by = "name") %>% 
  mutate(id = parse_integer(str_replace(Fold, "^Fold", "")))


# copy object to each cluster node
# cluster <- cluster %>% 
#   cluster_copy(tidyCV) %>% 
#   cluster_copy(df)

# partition data set to clusters
cvDF <- cvDF %>% 
  # partition(name, Fold, cluster = cluster) %>%
  group_by(name, Fold) %>% 
  # fit model on training part
  # fit model and save estimates in tidy format
  mutate(
    tidy_model = map2(Fold, design, .f = tidyer_fitter, 
                      tidyCV = tidyCV, data = df),
    pred = pmap(
      list(
        map(Fold, tidy_assessment, data = df, tidyCV = tidyCV),
        design,
        map(tidy_model, "estimate")
      ),
      sevenC:::predLogit
    ),
    label = map(map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), "loop")
  ) %>% 
  # collect results from cluster
  # collect() %>% 
  # ungroup
  ungroup()

write_rds(cvDF, path = paste0(outPrefix, "cvDF.rds"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF.rds"))

# gi <- read_rds(paste0(outPrefix, ".gi.rds"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF.rds"))
# designDF <- read_rds(paste0(outPrefix, "designDF.rds"))
#===============================================================================
# Performance Evaluation
#===============================================================================

# remove TF_only models
designDF <- designDF %>%
  filter(!str_detect(name, ".*_only$") ) %>%
  mutate(name = factor(name, name))

cvDF <- cvDF %>%
  filter(!str_detect(name, ".*_only$") )

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = cvDF$pred,
  labels = cvDF$label,
  modnames = as.character(cvDF$name),
  dsids = cvDF$id,
  posclass = levels(cvDF$label[[1]])[2],
  x_bins = 100)

write_rds(curves, paste0(outPrefix, ".curves.rds"))
# curves <- read_rds(paste0(outPrefix, ".curves.rds"))

# get data.frame with auc values
aucDF <-  as_tibble(auc(curves)) %>% 
  mutate(modnames = factor(modnames, designDF$name)) %>% 
  arrange(modnames)

aucDFmed <- aucDF %>%
  group_by(modnames, curvetypes) %>% 
  summarize(
    aucs_median = median(aucs, na.rm = TRUE),
    aucs_mean = mean(aucs, na.rm = TRUE),
    aucs_sd = sd(aucs, na.rm = TRUE)
  ) %>% 
  ungroup() 

write_feather(aucDFmed, paste0(outPrefix, ".aucDFmed.feather"))
# aucDFmed <- read_feather(paste0(outPrefix, ".aucDFmed.feather"))

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = modnames, y = aucs_mean, fill = modnames)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, vjust = 1.5) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.barplot.pdf"), w = 3.5, h = 7)


# Only AUC PRC as barplot
p <- ggplot(filter(aucDFmed, curvetypes == "PRC"), 
            aes(x = modnames, y = aucs_mean, fill = modnames)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance\n(auPRC)")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_PRC.barplot.pdf"), w = 3.5, h = 5)


#-------------------------------------------------------------------------------
# get ROC and PRC plots
#-------------------------------------------------------------------------------

subsetList = list(
  dist_orientation_motif = designDF$name[1:3],
  genomic = designDF$name[1:4],
  genomic_TF = designDF$name[1:10],
  all_TF = designDF$name[1:11],
  all = designDF$name
)

for (subStr in names(subsetList)){
  
  aucDFmedSub <- aucDFmed %>% 
    filter(modnames %in% subsetList[[subStr]])
  
  curveDF <- as_tibble(as.data.frame(curves)) %>% 
    filter(modname %in% subsetList[[subStr]]) %>% 
    mutate(modname = factor(modname, designDF$name))
  
  aucDFroc <- aucDFmedSub %>% 
    filter(curvetypes == "ROC")
  
  
  rocDF <- curveDF %>% 
    filter(type == "ROC")
  
  g <- ggplot(rocDF, aes(x = x, y = y, color = modname)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, lty = "dotted") +
    theme_bw() + theme(aspect.ratio = 1) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    scale_color_manual("",
                       values = designDF$color,
                       labels = paste0(
                         aucDFroc$modnames, 
                         ": AUC=", 
                         signif(aucDFroc$aucs_mean,3)
                       ),
                       guide = guide_legend(
                         override.aes = list(size = 2))
    ) +
    theme(text = element_text(size = 15), legend.position = c(.7,.4),
          legend.background = element_rect(fill = alpha('white', 0.1)))
  ggsave(g, file = paste0(outPrefix, ".", subStr, ".ROC.pdf"), w = 5, h = 5)
  
  #-------------------------------------------------------------------------------
  # get PRC plots
  #-------------------------------------------------------------------------------
  aucDFprc <- aucDFmedSub %>% 
    filter(curvetypes == "PRC")
  
  prc_base = precrec:::.get_pn_info(curves)$prc_base
  write_rds(prc_base, paste0(outPrefix, "prc_base.rds"))
            
  prcDF <- curveDF %>% 
    filter(type == "PRC")
  
  g <- ggplot(prcDF, aes(x = x, y = y, color = modname)) +
    geom_line() +
    geom_abline(intercept = prc_base, slope = 0, lty = "dotted") +
    theme_bw() + theme(aspect.ratio=1) +
    labs(x = "Recall", y = "Precision") +
    scale_color_manual(values = designDF$color) + 
    theme(text = element_text(size = 15), legend.position = "none") +
    ylim(0, 1)
  
  ggsave(g, file = paste0(outPrefix, ".", subStr, ".PRC.pdf"), w = 5, h = 5)
  

  # Only AUC PRC as barplot
  p <- ggplot(aucDFprc, 
              aes(x = modnames, y = aucs_mean, fill = modnames)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                  width = .25, position = position_dodge(width = 1)) + 
    geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
              size = 3, vjust = 1.5, angle = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), 
          legend.position = "none",
          text = element_text(size = 15)) +
    scale_fill_manual(values = designDF$color) +
    labs(x = "Models", y = "Prediction performance\n(AUC PRC)") +
    ylim(0, 0.45)
  # p
  ggsave(p, file = paste0(outPrefix, ".", subStr, ".AUC_ROC_PRC.barplot.pdf"), w = 3.5, h = 7)
  
}

#*******************************************************************************
# Binary prediction ------
#*******************************************************************************
# gi <- read_rds(paste0(outPrefix, ".gi.rds"))

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

# chr22 <- GRanges(seqinfo(gi))["chr22"]
# chr22GI <- subsetByOverlaps(gi, chr22, ignore.strand = TRUE)

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


# # write all pairs (with labeld true ones)
# writeLongRangeFormat(
#   gi = gi, 
#   score_vec = ifelse(mcols(gi)$loop == "Loop", 1, -1), 
#   output_file = paste0(outPrefix, ".gi.loop.longrange.txt")
)

#===============================================================================
#===============================================================================
