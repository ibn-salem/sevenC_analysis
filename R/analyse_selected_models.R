################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
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
require(ROCR)         # for binary clasification metrices

source("R/chromloop.functions.R")


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = parallel::detectCores() - 2

MIN_MOTIF_SIG <- 6
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

lfcPrefix <- file.path("results", paste0("v04_screen_TF_lfc.", 
                                         paste0("motifSig", MIN_MOTIF_SIG), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))


outPrefix <- file.path("results", paste0("v04_selected_models.", 
                                         paste0("motifSig", MIN_MOTIF_SIG), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

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
  cluster_library(packages = c("chromloop", "tidyverse"))

# evaluate help function code on each cluster
cluster_eval(cluster, source("R/chromloop.functions.R"))


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

# meta <- meta[1:3, ]

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  ancGR <- chromloop::motif.hg19.CTCF
  
  # filter for p-valu <= MIN_MOTIF_SIG
  ancGR <- ancGR[ancGR$sig >= MIN_MOTIF_SIG]
  
  seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
  
  # get all pairs within 1M distance
  gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)
  
  # add strand combinations
  gi <- chromloop::addStrandCombination(gi)
  
  # add motif score
  gi <- chromloop::addMotifScore(gi, colname = "sig")
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  trueLoopsTang2015 <- do.call(
    "c",
    lapply(LoopTang2015_GM12878_Files, 
           chromloop::parseLoopsTang2015, 
           seqinfo = seqInfo))
  
  gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
  gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))
  
}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Annotae with coverage and correlation -------------------------------------

if (!GI_LOCAL ) {
  
  # iterate over all ChIP-seq sata sets
  for (i in seq_len(nrow(meta))) {
    
    # check if sample is not present yet.
    if ( !paste0("cor_", meta$name[i]) %in% names(mcols(gi))){
      
      message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
      
      # add coverage  
      regions(gi) <- chromloop::addCovToGR(
        regions(gi), 
        meta$filePath[i], 
        window = WINDOW_SIZE,
        bin_size = BIN_SIZE,
        colname = paste0("cov_", meta$name[i])
      )
      
      # add correlations
      gi <- chromloop::applyToCloseGI(
        gi, 
        datcol = paste0("cov_", meta$name[i]),
        fun = cor, 
        colname = paste0("cor_", meta$name[i])
      )  
      
    }
  }  
  # Annotae with correlation across TFs -------------------------------------
  
  # get vector with coverage in whole anchor regions
  covCols <- paste0("cov_", meta$name)
  covDF <- as_tibble(as.data.frame(mcols(regions(gi))[,covCols])) %>% 
    dplyr::mutate_all(.funs = function(l) map_dbl(l, sum))
  
  mcols(regions(gi))[, "cov_sum"] <- NumericList(as_tibble(t(covDF)))
  
  gi <- chromloop::applyToCloseGI(
    gi, 
    datcol = "cov_sum",
    fun = cor, 
    colname = "across_TFs"
  )
  
  # save file for faster reload
  # save(gi, file = paste0(outPrefix, ".gi.Rdata"))
  write_rds(gi, paste0(outPrefix, ".gi.rds"))
  
} else {
  # load(paste0(outPrefix, ".gi.Rdata"))  
  gi <- read_rds(paste0(outPrefix, ".gi.rds"))
}


#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------

df <- as_tibble(as.data.frame(mcols(gi))) %>%
  mutate(
    id = 1:nrow(.),
    loop = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop"))
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
cluster <- cluster %>% 
  cluster_copy(tidyCV) %>% 
  cluster_copy(df)

# partition data set to clusters
cvDF <- cvDF %>% 
  partition(name, Fold, cluster = cluster) %>% 
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
      chromloop::pred_logit
    ),
    label = map(map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), "loop")
  )%>% 
  # collect results from cluster
  collect() %>% 
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
  modnames = cvDF$name,
  dsids = cvDF$id,
  posclass = levels(cvDF$label[[1]])[2],
  x_bins = 100)

write_rds(curves, paste0(outPrefix, ".curves.rds"))

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

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = modnames, y = aucs_mean, fill = modnames)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, hjust = 1, angle = 90, position = position_dodge(width = 1)) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
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
            size = 3, vjust = 1.5, angle = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance\n(AUC PRC)")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.barplot.pdf"), w = 3.5, h = 7)



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
    theme_bw() + theme(aspect.ratio=1) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    scale_color_manual("",
                       values = designDF$color,
                       labels = paste0(
                         aucDFroc$modnames, 
                         ": AUC=", 
                         signif(aucDFroc$aucs_mean,3)
                       ),
                       guide = guide_legend(
                         override.aes = list(size = 2),
                         reverse = TRUE)
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
    scale_color_manual(values = designDF$color,
                       labels = paste0(
                         aucDFroc$modnames, 
                         ": AUC=", 
                         signif(aucDFprc$aucs_mean,3)
                       ),
                       guide = guide_legend(
                         override.aes = list(size = 2),
                         reverse = TRUE)
    # ) + theme(legend.position = c(.75,.4))
    ) + theme(text = element_text(size = 15), legend.position = "none") +
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

  #-------------------------------------------------------------------------------
# Binary prediction
#-------------------------------------------------------------------------------



#===============================================================================
#===============================================================================
