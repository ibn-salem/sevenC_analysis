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

TRUE_LOOPS <- "HIC_ChIAPET"
# TRUE_LOOPS <- "HIC_ChIAPET_CaptureC"

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
COL_TF <- grDevices::colors()[str_detect(grDevices::colors(), "^((?!(gr(a|e)y|white)).)*[^\\d]$")]
# COL_TF <- COL_TF[!COL_TF %in% c("white", "black", "aliceblue", "azure", "beige", "bisque", "cornsilk", "cyan", "darkorchid", "coral", "darkmagenta")]
# pie(rep(1, length(COL_TF)), col = COL_TF, labels = COL_TF, main = length(COL_TF))
# COL_TF = colorRamps::primary.colors(124)

COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)]
names(COL_LOOP) <- c("No loop", "Loop")

outPrefix <- file.path("results", paste0("v04_screen_TF_lfc.", 
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
# metaFile <- "data/ENCODE/metadata.fltBam.tsv"
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
  mutate(name = paste0(TF, "_lfc")) %>%
  select(TF, name, filePath, everything())

# adapte colors
COL_TF <- COL_TF[1:nrow(meta)]
names(COL_TF) <- meta$TF
COL_TF[match(SELECTED_TF, names(COL_TF))] <- COL_SELECTED_TF_2
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
  

  # ol <- IRanges::overlapsAny(gi, trueLoops)
  # gi$Loop_Rao_GM12878 <- factor(ol, c(FALSE, TRUE), c("No loop", "Loop"))
  
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
    colname = "cor_across_TFs"
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

# save(df, file = paste0(outPrefix, ".df.Rdata")) 
write_feather(df, paste0(outPrefix, ".df.feather"))
# load(paste0(outPrefix, ".df.Rdata"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

#-------------------------------------------------------------------------------
# Compare dist vs. log10(dist) in prediction 
#-------------------------------------------------------------------------------
formula_list <- list(
  dist = as.formula("loop ~ dist + strandOrientation + score_min + cor_CTCF_lfc"),
  log10dist = as.formula("loop ~ dist_log10 + strandOrientation + score_min + cor_CTCF_lfc")
)

mod_list <- formula_list %>% 
  map(glm, family = binomial(), data = df)

param_list <- mod_list %>% 
  map(tidy) %>% 
  map("estimate")

pred_list <- map2(formula_list, param_list, 
                  ~ chromloop::pred_logit(data = df, formula = .x, betas = .y))


mdat <- mmdata(pred_list, labels = list(df$loop, df$loop), modnames = c("dist", "log10_dist"))
curves <- evalmod(mdat)

p <- autoplot(curves)
ggsave(p, paste0(outPrefix, ".compare_dist_vs_log10dist.corCTCF_lfc.pdf"))
#-------------------------------------------------------------------------------
# split whole data set into EDA and prediction 
#-------------------------------------------------------------------------------

# randomly sample 10% of the data for EDA analysis
# split_eda <- sample.int(nrow(df), 0.1 * nrow(df))
# 
# edaDF <- df[split_eda,]
# modDF <- df[-split_eda,]
# 
# write_feather(edaDF, paste0(outPrefix, ".edaDF.feather"))
# write_feather(modDF, paste0(outPrefix, ".modDF.feather"))
# 
# #remove all but 3 TF columns
# rmNames <- paste0("cor_", meta$name[4:nrow(meta)])
# df <- df %>%
#   select(-match(rmNames, names(.)))

#===============================================================================
# Training and Cross-validation
#===============================================================================

# useTF <- meta$name
#filter for a subset of TFS
#

# k fold cross validation with 1 repeats
set.seed(3579)

tidyCV <- df %>% 
  vfold_cv(V = K, repeats = 1) %>% 
  tidy()

write_feather(tidyCV, paste0(outPrefix, ".tidyCV.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))

# get design formula for each TF
designDF <- tibble(
  name = meta$name,
  design = map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x)) )
)

# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = meta$name, Fold) %>% 
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
                      tidyCV = tidyCV, data = df)
  )%>% 
  # collect results from cluster
  collect()

write_rds(cvDF, path = paste0(outPrefix, "cvDF_trained.rds"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF_trained.rds"))

# df <- read_feather(paste0(outPrefix, ".df.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF_trained.rds"))

#===============================================================================
# Prediction
#===============================================================================

# add prediction using individual TF specific models
cvDF <- cvDF %>%
  partition(name, Fold, cluster = cluster) %>% 
  mutate(
    pred_specificTF = pmap(
      list(
        map(Fold, tidy_assessment, data = df, tidyCV = tidyCV),
        design,
        map(tidy_model, "estimate")
      ),
      chromloop::pred_logit
    ),
    label = map(map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), "loop")
  ) %>% 
  collect()

# ungroup and get TF as seprate column  
cvDF <- cvDF %>% 
  ungroup() %>% 
  mutate(TF = str_replace(name, "_lfc", ""))

# save with predictions
write_rds(cvDF, path = paste0(outPrefix, "cvDF_specificTFpred.rds"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF_specificTFpred.rds"))

#===============================================================================
# Performance Evaluation
#===============================================================================

# gather the two differnt prediction types into one colum
# extract only the needed columns
evalDF <- cvDF %>% 
  ungroup() %>% 
  mutate(
    fold = id,
    pred = pred_specificTF,
    modnames = TF
  ) %>% 
  select(modnames, fold, pred, label)

# save evalDF
# write_rds(evalDF, paste0(outPrefix, ".evalDF_specificTF.rds"))  
# evalDF <- read_rds(paste0(outPrefix, ".evalDF_specificTF.rds"))

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = evalDF$pred,
  labels = evalDF$label,
  modnames = evalDF$modnames,
  dsids = evalDF$fold,
  posclass = levels(evalDF$label[[1]])[2],
  x_bins = 100)


# get data.frame with auc values
aucDF <-  as_tibble(auc(curves))

# get ranked modle names
ranked_models <- aucDF %>% 
  filter(curvetypes == "PRC") %>% 
  group_by(modnames) %>% 
  summarize(
    auc_mean = mean(aucs, na.rm = TRUE)
  ) %>% 
  arrange(desc(auc_mean)) %>% 
  pull(modnames)

write_rds(ranked_models, paste0(outPrefix, "ranked_models.rds"))

# order aucDF by ranks
aucDF <- aucDF %>% 
  mutate(modnames = factor(modnames, ranked_models)) %>% 
  arrange(modnames)

write_feather(aucDF, paste0(outPrefix, "aucDF_specificTF.feather"))
# aucDF <- read_feather(paste0(outPrefix, "aucDF_specificTF.feather"))
# ranked_models <- levels(aucDF$modnames)[1:N_TOP_MODELS]

#-------------------------------------------------------------------------------
# Take parameters from best N models
#-------------------------------------------------------------------------------
TFspecific_ModelDF <- cvDF %>%
  ungroup() %>% 
  mutate(TF = str_replace(name, "_lfc", "")) %>% 
  select(TF, id, tidy_model) %>%
  unnest(tidy_model) %>% 
  # rename cor_* terms to only cor
  mutate(term = str_replace(term, "^cor_.*", "cor")) %>% 
  mutate(term = factor(term, unique(term))) %>% 
  group_by(TF, term) %>% 
  summarize(
    estimate_mean = mean(estimate, na.rm = TRUE),
    estimate_median = median(estimate, na.rm = TRUE),
    estimate_sd = sd(estimate, na.rm = TRUE)
  ) %>% 
  mutate(term = parse_character(term))

write_tsv(TFspecific_ModelDF, paste0(outPrefix, ".TFspecific_ModelDF.tsv"))



# take mean/meidan of estimates across TFs and folds
allTfModelDF <- cvDF %>% 
  select(TF, id, tidy_model) %>%
  unnest(tidy_model) %>% 
  # rename cor_* terms to only cor
  mutate(term = str_replace(term, "^cor_.*", "cor")) %>% 
  mutate(term = factor(term, unique(term))) %>% 
  group_by(term) %>% 
  summarize(
    estimate_mean = mean(estimate, na.rm = TRUE),
    estimate_median = median(estimate, na.rm = TRUE),
    estimate_sd = sd(estimate, na.rm = TRUE)
  ) %>% 
  mutate(term = parse_character(term))

write_tsv(allTfModelDF, paste0(outPrefix, ".allTfModelDF.tsv"))

# combine models to a single one
bestNModelDF <- cvDF %>%
  filter(TF %in% ranked_models[1:N_TOP_MODELS]) %>%
  select(TF, id, tidy_model) %>% 
  unnest(tidy_model) %>% 
  # rename cor_* terms to only cor
  mutate(term = str_replace(term, "^cor_.*", "cor")) %>% 
  mutate(term = factor(term, unique(term))) %>% 
  group_by(term) %>% 
  summarize(
    estimate_mean = mean(estimate, na.rm = TRUE),
    estimate_median = median(estimate, na.rm = TRUE),
    estimate_sd = sd(estimate, na.rm = TRUE)
  ) %>% 
  mutate(term = parse_character(term))

write_tsv(bestNModelDF, paste0(outPrefix, ".bestNModelDF.tsv"))
# bestNModelDF <- read_tsv(paste0(outPrefix, ".bestNModelDF.tsv"))

#-------------------------------------------------------------------------------
# Plot parameters
#-------------------------------------------------------------------------------

# get tidy model DF
modelDF <- cvDF %>% 
  unnest(tidy_model, .drop = TRUE) %>% 
  mutate(term = str_replace(term, "^cor_.*", "cor")) %>% 
  mutate(term = factor(term, unique(term))) %>% 
  filter(TF %in% SELECTED_TF) %>% 
  select(name, id, term, estimate)

allTF <- allTfModelDF %>%
  mutate(estimate = estimate_mean, 
         name = "allTF") %>% 
  select(name, term, estimate)

bestN <- bestNModelDF %>% 
  mutate(estimate = estimate_mean,
         name = "bestN") %>% 
  select(name, term, estimate)

modelDF <- modelDF %>% 
  bind_rows(allTF, bestN)

write_tsv(modelDF, paste0(outPrefix, ".modelDF.tsv"))

# # add meta data
# modelDF <- modelDF %>% 
#   left_join(meta, by = "name")

paramByModel <- modelDF %>% 
  group_by(name, term) %>% 
  summarize(
    n = n(),
    estimate_mean = mean(estimate),
    estimate_sd = sd(estimate)
  )

p <- ggplot(paramByModel, aes(x = name, y = estimate_mean, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(
    ymin = estimate_mean - estimate_sd, 
    ymax = estimate_mean + estimate_sd), width = 0.25) +
  geom_text(aes(label = round(estimate_mean, 2)), hjust = "inward") + 
  facet_grid(. ~ term , scales = "free_x") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + #scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, file = paste0(outPrefix, ".selected_models.paramter.barplot.pdf"), w = 12, h = 6)


#-------------------------------------------------------------------------------
# Add predictions
#-------------------------------------------------------------------------------

# Prediction using custom model

Rad21mod <- cvDF %>%
  filter(TF == "RAD21", id == 1) %>%
  pull(tidy_model)

Rad21mod <- Rad21mod[[1]]

# defaultDesign <- as.formula("loop ~ dist + strandOrientation + score_min + cor")

# copy object to each cluster node
cluster <- cluster %>% 
  cluster_copy(allTfModelDF) %>% 
  cluster_copy(bestNModelDF) %>% 
  cluster_copy(Rad21mod)
  
cvDF <- cvDF %>% 
  partition(name, Fold, cluster = cluster) %>%
  mutate(
    # add prediction using single TF model
    pred_allTF = map2(
      .x = map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), 
      .y = design, 
      .f = pred_logit,
      betas = allTfModelDF$estimate_mean),
    # add prediction using N best models
    pred_bestNTF = map2(
      .x = map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), 
      .y = design, 
      .f = pred_logit,
      betas = bestNModelDF$estimate_mean),
    pred_rad21 = map2(
      .x = map(Fold, tidy_assessment, data = df, tidyCV = tidyCV),
      .y = design,
      .f = pred_logit,
      betas = Rad21mod$estimate)
    ) %>% 
  collect()

# save with predictions of n best models
write_rds(cvDF, path = paste0(outPrefix, "cvDF_withPred.rds"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF_withPred.rds"))

#-------------------------------------------------------------------------------
# AUC of different predictions
#-------------------------------------------------------------------------------

# gather the two differnt prediction types into one colum
evalDF <- cvDF %>% 
  ungroup() %>% 
  gather(starts_with("pred_"), key = "pred_type", value = "pred") %>% 
  mutate(
    fold = id,
    modnames = paste0(TF, "_", str_replace(pred_type, "pred_", ""))
  ) %>% 
  select(modnames, fold, pred, label) 

# save evalDF
write_rds(evalDF, paste0(outPrefix, ".evalDF.rds"))  
# evalDF <- read_rds(paste0(outPrefix, ".evalDF.rds"))

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = evalDF$pred,
  labels = evalDF$label,
  modnames = evalDF$modnames,
  dsids = evalDF$fold,
  posclass = levels(evalDF$label[[1]])[2],
  x_bins = 100)


# get data.frame with auc values
aucDF <- as_tibble(auc(curves)) %>% 
  separate(modnames, into = c("name", "pred_type"), sep = "_", remove = FALSE) %>% 
  mutate(name = factor(name, ranked_models))

write_feather(aucDF, paste0(outPrefix, "aucDF.feather"))
# aucDF <- read_feather(paste0(outPrefix, "aucDF.feather"))

# build color vector with TFs as names
# COL_TF <- COL_TF[seq(1, length(ranked_models))]
# names(COL_TF) <- ranked_models

aucDFmed <- aucDF %>%
  group_by(name, pred_type, curvetypes) %>% 
  summarize(
    aucs_median = median(aucs, na.rm = TRUE),
    aucs_mean = mean(aucs, na.rm = TRUE),
    aucs_sd = sd(aucs, na.rm = TRUE)
  ) %>% 
  mutate(topN = factor(name %in% ranked_models[1:N_TOP_MODELS], c(TRUE, FALSE), c("top_n", "rest")))

# aucDFmed <- aucDFmed %>%
#   filter(name %in% ranked_models[1:6])

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = name, y = aucs_mean, fill = pred_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, hjust = 1, angle = 90, position = position_dodge(width = 1)) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom") +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Models", y = "AUC")

ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF_and_predType.barplot.pdf"), w = 14, h = 7)

# plot only selected TFs
selectedDFmed <- aucDFmed %>% 
  filter(name %in% SELECTED_TF)

p <- ggplot(selectedDFmed, aes(x = name, y = aucs_mean, fill =  pred_type)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, hjust = 1, angle = 90, position = position_dodge(width = 1)) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom") +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Models", y = "AUC")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF_and_predType_selectedTF.barplot.pdf"), w = 7, h = 7)

# boxplot of AUCs across TFs
p <- ggplot(aucDFmed, aes(x = pred_type, y = aucs_mean, fill = pred_type)) +
  geom_boxplot() +
  facet_grid(curvetypes ~ topN, margins = "topN") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom") +
  scale_fill_brewer(palette = "Blues") + 
  labs(x = "Models", y = "AUC")

ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF_and_predType.boxplot.pdf"), w = 7, h = 7)


#-------------------------------------------------------------------------------
# Barplot of AUC for only specific model
#-------------------------------------------------------------------------------
# aucDF <- read_feather(paste0(outPrefix, "aucDF.feather"))

# build color vector with TFs as names
# COL_TF <- COL_TF[seq(1, length(ranked_models))]
# names(COL_TF) <- ranked_models

aucDFcombined <- aucDF %>%
  # take only the prediction value of the TF specific model
  filter(pred_type == "specificTF") %>% 
  group_by(name, curvetypes) %>% 
  summarize(
    aucs_median = median(aucs, na.rm = TRUE),
    aucs_mean = mean(aucs, na.rm = TRUE),
    aucs_sd = sd(aucs, na.rm = TRUE)
  )
  
#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFcombined, aes(x = name, y = aucs_mean, 
                          ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd, 
                          fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = .25) + 
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, hjust = 1, angle = 90, position = position_dodge(width = 1)) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "TF", y = "AUC")

ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)

# only PRC AUC
prcDF <- aucDFcombined %>% 
  filter(curvetypes == "PRC")

p <- ggplot(prcDF, aes(x = name, y = aucs_mean, 
                               ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd, 
                               fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = .25) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "", y = "Prediction perfromance\n(AUC PRC)")

ggsave(p, file = paste0(outPrefix, ".AUC_PRC.by_TF.barplot.pdf"), w = 14, h = 3.5)


#===============================================================================
# TODO: Did run it only until here.
#===============================================================================
stop("Stop here!")

#-------------------------------------------------------------------------------
# get ROC plots
#-------------------------------------------------------------------------------

# get AUC of ROC and PRC curves for all 
evalDFsub <- evalDF %>% 
  mutate(TF = str_split_fixed(modnames, "_", 2)[ ,1]) %>%
  mutate(type = str_split_fixed(modnames, "_", 2)[ ,2]) %>%
  filter(TF %in% SELECTED_TF) %>% 
  filter(type == "specificTF")

curvesSub <- evalmod(
  scores = evalDFsub$pred,
  labels = evalDFsub$label,
  modnames = evalDFsub$modnames,
  dsids = evalDFsub$fold,
  posclass = levels(evalDFsub$label[[1]])[2],
  x_bins = 100)


aucDFroc <- auc(curvesSub) %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curvesSub, "ROC", show_cb = TRUE) + 
  scale_color_manual(values = COL_SELECTED_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFroc$aucs,3)
                     ),
                     guide = guide_legend(
                       override.aes = list(size = 2),
                       reverse = TRUE)
  ) + 
  theme(legend.position = c(.75,.4))

# g
ggsave(g, file= paste0(outPrefix, ".selectedTF_specificTF.ROC.pdf"), w = 5, h = 5)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(curves, "PRC", size = 4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFprc$modnames, 
                       ": AUC=", 
                       signif(aucDFprc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position = "none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)

#-------------------------------------------------------------------------------
# Binary prediction
#-------------------------------------------------------------------------------
require("ROCR")

evalDF <- cvDF %>% 
  select(name, id, pred_allTF, label) %>%
  # filter(name %in% SELECTED_TF) %>% 
  filter(id == 1)


# get prediction object with measurements
predObj <- ROCR::prediction(evalDF$pred_allTF, evalDF$label, levels(evalDF$label[[1]]))
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
  name = rep(evalDF$name, times = map_int(slot(perfObj, "x.values"), length))
)

p <- ggplot(f1DF, aes(x = cutoff, y = f1_score, color = name)) +
  geom_line() +
  theme_bw() + theme(legend.position = "none") 
ggsave(p, file = paste0(outPrefix, ".allTF.pred_allTF.f1-score_vs_cutoff.pdf"), w = 5, h = 5)

# plot curve only for selected TFs
p <- ggplot(filter(f1DF, name %in% SELECTED_TF), aes(x = cutoff, y = f1_score, color = name)) +
  geom_line() +
  theme_bw() + theme(legend.position = "bottom") + 
  scale_color_manual(values = COL_SELECTED_TF)
ggsave(p, file = paste0(outPrefix, ".selectedTF.pred_allTF.f1-score_vs_cutoff.pdf"), w = 5, h = 5)

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
