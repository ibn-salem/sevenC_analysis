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
library(ROCR)         # for binary clasification metrices
library(furrr)        # for parallelization

plan(multiprocess, workers = 3)
plan()

source("R/sevenC.functions.R")

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(10, parallel::detectCores() - 1)

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
    message("INFO: File: ", meta$filePath[[i]], " Exists: ", file.exists(meta$filePath[[i]]))
    
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
  
  # compute sum of signal at first and second anchor as separate predictor colums
  
  reg <- mcols(regions(gi))
  idx1 <- anchorIds(gi, "first")
  idx2 <- anchorIds(gi, "second")
  
  for (TF in meta$name) {
    
    message("INFO: --> Working on sample: ", TF, " <--")
    
    # get sum of all position values as total signal per anchor region
    reg[, str_c("sum_", TF)] <- reg[, TF] %>% as.list() %>% map_dbl(sum)
    
    # add first and second anchor signal separatly to gi object
    mcols(gi)[, str_c("cov1_", TF)] <-  reg[, str_c("sum_", TF)][idx1]
    mcols(gi)[, str_c("cov2_", TF)] <-  reg[, str_c("sum_", TF)][idx2]
  }
  
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
    paste0(meta$name, "_cov"),
    meta$name
  ),
  design = c(
    designList,
    map(meta$name, ~as.formula(paste0("loop ~ cor_", .x)) ),
    map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cov1_", .x, " + cov2_", .x)) ),
    map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x)) )
  ),
  color = c(
    c("greenyellow", "gold2", "khaki", "brown4", "darkgray", "grey30"),
    COL_SELECTED_TF_1,
    COL_SELECTED_TF_1,
    COL_SELECTED_TF_2
  )
) %>% 
  mutate(name = factor(name, c(names(designList)[1:4], paste0(meta$name, "_only"), paste0(meta$name, "_cov"), meta$name, "all_TF", "across_TFs"))) %>% 
  arrange(name)

write_rds(designDF, paste0(outPrefix, "designDF.rds"))

# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = designDF$name, Fold) %>% 
  # add design formular for each TF
  left_join(designDF, by = "name") %>% 
  mutate(id = parse_integer(str_replace(Fold, "^Fold", "")))

# partition data set to clusters
cvDF <- cvDF %>% 
  group_by(name, Fold) %>% 
  # fit model on training part
  # fit model and save estimates in tidy format
  mutate(
    tidy_model = future_map2(Fold, design, .f = tidyer_fitter, 
                      tidyCV = tidyCV, data = df),
    pred = future_map(
      list(
        map(Fold, tidy_assessment, data = df, tidyCV = tidyCV),
        design,
        map(tidy_model, "estimate")
      ),
      sevenC:::predLogit
    ),
    label = map(map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), "loop")
  ) %>% 
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
  ungroup() %>% 
  mutate(
    curvetypes = factor(curvetypes, c("ROC", "PRC"))
  )

write_feather(aucDFmed, paste0(outPrefix, ".aucDFmed.feather"))
# aucDFmed <- read_feather(paste0(outPrefix, ".aucDFmed.feather"))

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = modnames, y = aucs_mean, fill = modnames)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) +
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, vjust = 1.5) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.barplot.pdf"), w = 5, h = 7)


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

for (subStr in names(subsetList)) {
  
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
    theme_bw() + theme(aspect.ratio = 1) +
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
  
  ggsave(p, file = paste0(outPrefix, ".", subStr, ".AUC_ROC_PRC.barplot.pdf"), w = 3.5, h = 7)
  
}

#===============================================================================
#===============================================================================
