################################################################################
# Predict chromatin loopsin HeLa cells and validate with Hi-C and ChIA-PET
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(rtracklayer)  # to import() BED files
require(rsample)
require(pryr) # for object_size()
require(feather)      # for efficient storing of data.frames
require(multidplyr)   # for partition() and collect() to work in parallel

source("R/chromloop.functions.R")


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

# meta <- meta[1:3, ]

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
                           ~ chromloop:::predLogit(df, formula = .x, betas = .y)),
    pred_allTF = map(design, ~ chromloop:::predLogit(df, .x, allTfModelDF$estimate_mean)),
    pred_bestN = map(design, ~ chromloop:::predLogit(df, .x, bestNModelDF$estimate_mean))
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
p <- ggplot(aucDF, aes(x = name, y = aucs, fill = name)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = round(aucs, 2), y = aucs), size = 3, vjust = 1.2, angle = 0) +
  facet_grid(curvetypes ~ pred_type, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_SELECTED_TF_2) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC_byPredtype.barplot.pdf"), w = 7, h = 7)


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
  theme_bw() + theme(aspect.ratio=1) +
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
  theme(text = element_text(size = 15), legend.position = c(.7,.3),
        legend.background = element_rect(fill = alpha('white', 0.1)))
ggsave(g, file = paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

#-------------------------------------------------------------------------------
# get PRC plots
#-------------------------------------------------------------------------------
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC", pred_type == "specificTF")

# prc_base = precrec:::.get_pn_info(curves)$prc_base
# write_rds(prc_base, paste0(outPrefix, "prc_base.rds"))

prcDF <- curveDFsub %>% 
  filter(type == "PRC", pred_type == "specificTF")

g <- ggplot(prcDF, aes(x = x, y = y, color = name)) +
  geom_line() +
  geom_abline(intercept = prc_base, slope = 0, lty = "dotted") +
  theme_bw() + theme(aspect.ratio = 1) +
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
  theme(legend.position = "none")

ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)


#-----------------------------------------------------------------------------
# Binary prediction
#-----------------------------------------------------------------------------
stop("Stop here. Rest needs to be implemented")

specificTFeval <- evalDF %>% 
  filter(pred_type == "pred_specificTF") %>% 
  mutate(label = list(df$loop))

# get prediction object with measurements
predObj <- ROCR::prediction(specificTFeval$pred, specificTFeval$label, levels(evalDF$label[[1]]))
perfObj <- ROCR::performance(predObj, measure = "f")

# get cutoff and f1-score for each model as lists
cutoff_List <- slot(perfObj, "x.values")
f1_List <- slot(perfObj, "y.values")

# # plot f1-score vs. cutoffs
# pdf(paste0(outPrefix, ".selectedTF.pred_allTF.f1-score_vs_cutoff.pdf"))
#   plot(performance(predObj, measure = "f"), col = brewer.pal(nrow(evalDF), "Dark2"))
# dev.off()

# plot f1-score vs. cutoffs using ggplot2
f1DF <- tibble(
  cutoff = unlist(cutoff_List),
  f1_score = unlist(f1_List),
  name = rep(unique(evalDF$name), map_int(cutoff_List, length))
)

p <- ggplot(f1DF, aes(x = cutoff, y = f1_score, color = name)) +
  geom_line() +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = COL_SELECTED_TF)
ggsave(p, file = paste0(outPrefix, ".pred_specificTF.f1-score_vs_cutoff.pdf"), w = 5, h = 5)

stop("Stop here. Rest needs to be implemented")

#-------------------------------------------------------------------------------
# get cutoff with maximal f1-score
#-------------------------------------------------------------------------------

f1ModelDF <- evalDF %>% 
  mutate(
    cutoffs = cutoff_List,
    f1_score = f1_List,
    max_idx = map_int(f1_score, which.max),
    max_cutoff = map2_dbl(cutoffs, max_idx, ~ .x[[.y]]),
    max_f1 = map2_dbl(f1_score, max_idx, ~ .x[[.y]])
  ) %>% 
  select(name, max_idx, max_cutoff, max_f1)

write_rds(f1ModelDF, paste0(outPrefix, ".f1ModelDF.rds"))  
write_tsv(f1ModelDF, paste0(outPrefix, ".f1ModelDF.tsv"))  

allTFf1ModelDF <- f1ModelDF %>% 
  summarize(
    mean_max_cutoff = mean(max_cutoff, na.rm = TRUE),
    median_max_cutoff = median(max_cutoff, na.rm = TRUE),
    sd_max_cutoff = sd(max_cutoff, na.rm = TRUE)
  )
write_tsv(allTFf1ModelDF, paste0(outPrefix, ".f1ModelDF.tsv"))  

# output the mean of the top N models
ranked_models <- levels(aucDF$modnames)[1:10]

topNf1ModelDF <- f1ModelDF %>% 
  filter(name %in% ranked_models) %>% 
  summarize(
    mean_max_cutoff = mean(max_cutoff, na.rm = TRUE),
    median_max_cutoff = median(max_cutoff, na.rm = TRUE),
    sd_max_cutoff = sd(max_cutoff, na.rm = TRUE)
  )

write_tsv(topNf1ModelDF, paste0(outPrefix, ".topNf1ModelDF.tsv"))  


#===============================================================================
#===============================================================================
