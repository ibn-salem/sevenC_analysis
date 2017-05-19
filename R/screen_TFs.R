################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(forcats)      # for factors
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(rtracklayer)  # to import() BED files
require(colorRamps)  # to pick different colors
require(scales)     # for scales in plotting and percent() formatin function
require(stringr)    # for string function and regular expressions
require(BiocParallel) # for parallelisation
require(biobroom)     # to make BioC classes tidy
require(grid)         # for textGrob
require(magrittr)     # for "tee" pipe %T>% 

#-----------------------------------------------------------------------
# Options for parallel computation
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed = 34312)
# set options
register(multicorParam)
# bpparam() # to print current options


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
MIN_MOTIF_SIG <- 6


WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation

TRUE_LOOPS <- "HIC_ChIAPET"
# TRUE_LOOPS <- "HIC_ChIAPET_CaptureC"

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
# pie(rep(1, length(COL_TF)), col=COL_TF, labels=COL_TF, main=length(COL_TF))

# COL_TF = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# str_view(grDevices::colors(), "[^\d]$")
# str_view(grDevices::colors(), "^((?!gr(a|e)y).)*[^\\d]$")
COL_TF <- grDevices::colors()[str_detect(grDevices::colors(), "^((?!(gr(a|e)y|white)).)*[^\\d]$")]
COL_TF <- COL_TF[!COL_TF %in% c("white", "black", "aliceblue", "azure", "beige", "bisque", "cornsilk", "cyan", "darkorchid", "coral", "darkmagenta")]
# pie(rep(1, 50), col=COL_TF, labels = COL_TF)
# pie(rep(1, length(COL_TF)), col=COL_TF)
# COL_TF = colorRamps::primary.colors(12)

#barplot(1:10, col=COL_TF)
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path(
  "results", 
  paste0("v02_screen_TF.",
           paste0("motifSig", MIN_MOTIF_SIG),
         "_w", WINDOW_SIZE,
         "_b", BIN_SIZE
  )
)

dir.create(dirname(outPrefix), showWarnings = FALSE)

# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_GM12878_Files <- c(
  "data/Tang2015/GSM1872886_GM12878_CTCF_PET_clusters.txt",
  "data/Tang2015/GSM1872887_GM12878_RNAPII_PET_clusters.txt")

CaptureHiC_Files <- c(
  "data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt",
  "data/Mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt"
)

metaFile <- "data/ENCODE/metadata.fcDF.tsv"

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

meta <- meta %>% 
  mutate(outType = str_replace(`Output type`, "fold change over control", "fold_change")) %>% 
  unite(TF, outType, sep = "_", col = "name", remove = FALSE) %>% 
  select(TF, outType, everything())


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
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  trueLoopsTang2015 <- do.call(
    "c",
    lapply(LoopTang2015_GM12878_Files, 
           chromloop::parseLoopsTang2015, 
           seqinfo = seqInfo))
  
  trueCaptureHiC <- do.call(
    "c",
    lapply(CaptureHiC_Files,
           chromloop::parseCaptureHiC, seqinfo = seqInfo)
  )
  trueCaptureHiC <- trueCaptureHiC[trueCaptureHiC$log_observed_expected >= 10]
  # trueCaptureHiC <- chromloop::parseCaptureHiC(
  #   CaptureHiC_Files[1], seqinfo = seqInfo)
  
  # ol <- IRanges::overlapsAny(gi, trueLoops)
  # gi$Loop_Rao_GM12878 <- factor(ol, c(FALSE, TRUE), c("No loop", "Loop"))
  
  gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
  gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
  gi <- addInteractionSupport(gi, trueCaptureHiC, "Loop_Mifsud2015_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))
  
}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Annotae with coverage and correlation -------------------------------------

if (!GI_LOCAL ) {
  
  for (i in seq_len(nrow(meta))) {
    
    message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
    
    stopifnot(file.exists(meta$filePath[i]))
    
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
  
  # Annotae with correlation across TFs -------------------------------------
  
  covCols <- paste0("cov_", meta$name)
  covDF <- as_tibble(as.data.frame(mcols(regions(gi))[,covCols])) %>% 
    mutate_all(.funs = purrr::map_dbl, .f = sum)
  
  mcols(regions(gi))[, "cov_sum"] <- NumericList(as_tibble(t(covDF)))
  
  gi <- chromloop::applyToCloseGI(
    gi, 
    datcol = "cov_sum",
    fun = cor, 
    colname = "cor_across_TFs"
  )
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.Rdata"))
  
} else {
  load(paste0(outPrefix, ".gi.Rdata"))  
}


#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------


df <- as_tibble(as.data.frame(mcols(gi))) %>%
  mutate(
    id = row_number(.),
    HIC_ChIAPET_CaptureC = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop" | Loop_Mifsud2015_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop")),
    HIC_ChIAPET = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop"))
  ) 


if (TRUE_LOOPS == "HIC_ChIAPET_CaptureC") {
  df <- df %>% mutate(loop =  HIC_ChIAPET_CaptureC)
  outPrefix <- str_c(outPrefix, "_HIC_ChIAPET_CaptureC")
} else {
  df <- df %>% mutate(loop =  HIC_ChIAPET)
}

df <- df %>% 
  select(id, loop, everything()) %>% 
  select(-matches("cor_across_TFs")) 

save(df, file = paste0(outPrefix, ".df.Rdata"))  
# load(paste0(outPrefix, ".df.Rdata"))


# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = name, value = cor) 

tidyDF <- tidyDF %>% 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(name = str_sub(name, 5)) %>% 
  separate(name, into = c("TF", "outType"), 
           sep = "_", extra = "merge", remove = FALSE)


# ------------------------------------------------------------------------------
# NAs in correlation ------------------------------------------------------
# ------------------------------------------------------------------------------

naDF <- tidyDF %>% 
  group_by(TF, outType, loop) %>% 
  summarise(
    n = n(),
    nNotNA = sum(!is.na(cor)),
    percentNotNA = nNotNA / n * 100
  )

write_tsv(naDF, paste0(outPrefix, ".naDF.tsv"))

p <- ggplot(naDF, aes(x = TF, y = percentNotNA, fill = TF)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = round(percentNotNA, 2)), vjust = 1.5) +
  facet_grid(loop + outType ~ ., scales = "free_x") + 
  theme_bw() + theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none") + 
  scale_fill_manual(values = COL_TF)
# p
ggsave(p, file = paste0(outPrefix, ".percentNotNA.by_TF_and_loop.barplot.pdf"), w = 14, h = 7)

# percent of non-missing observations
completeDF <- df %>% 
  select(-starts_with("Loop_")) %>% 
  mutate(complete = complete.cases(.)) %>% 
  count(complete)

# ------------------------------------------------------------------------------
# Compare Correlation with Boxplot
# ------------------------------------------------------------------------------

tidySubDF <- tidyDF %>% 
  # filter(TF %in% c("CTCF", "Stat1")) %>% 
  sample_n(min(nrow(tidyDF), 10^7))

p <- ggplot(tidySubDF, aes(x = TF, y = cor, col = loop)) +
  # geom_violin(aes(fill = TF), lwd = 1.5) + 
  geom_boxplot() +
  facet_grid(outType ~ ., scales = "free_x") + 
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    # text = element_text(size=20), 
    # axis.text.x=element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Correlation of ChIP-seq signal", x="")
# geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
# p
ggsave(p, file = paste0(outPrefix, ".cor.by_TF_and_loop.boxplot.pdf"), w = 14, h = 7)

# ------------------------------------------------------------------------------
# Predict loops --------------------------------------------------------
# ------------------------------------------------------------------------------

useTF <- meta$name

# dived in training and test set by taking 9/10 for training and 1/10 for testing
# testCases <- sample.int(nrow(df), size = round(nrow(df)/10))
# testCasesVec <- ifelse(1:nrow(df) %in% testCases, "train", "test")
# train <- df[-testCases,]
# test <- df[testCases,]

# split dataset in randomly in K equal sized parts
fold <- sample(cut(seq(1, nrow(df)), breaks = K, labels = FALSE))

# tidy prediction --------------------------------------------------
byTF <- tidyDF %>%
  mutate(fold = rep(fold, length(useTF))) %>% 
  group_by(name) %>% 
  nest()

singleTF_model <- function(df, fold_idx = 1) {
  
  glm(
    loop ~ dist + strandOrientation + cor, 
    family = binomial(link = 'logit'), 
    data = subset(df, fold != fold_idx)
    )  
}

add_predict_fold <- function(df, model, fold_idx = 1) {
  df[df$fold == fold_idx, "pred"] <- predict(
    model, 
    newdata = subset(df, fold == fold_idx), 
    type = "response")
  return(df)
}


# m <- singleTF_model(byTF$data[[1]], 2)
# models <- map(byTF$data, singleTF_model)

# for (k in 1:K) {
# }  

# fit model for each TF and get predictions
byTF <- byTF %>%
  mutate(model = map(byTF$data, singleTF_model, fold_idx = 1)) %>%
  mutate(data = map2(data, model, add_predict_fold, fold_idx = 1))

#-------------------------------------------------------------------------------
# get predictions from single TF with dist
#-------------------------------------------------------------------------------
# modelList <- lapply(useTF, function(name) {
#   
#   message("INFO: Fit model for TF: ", name)
#   
#   design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + cor_", name))
#   
#   model <- glm(
#     design, 
#     family = binomial(link = 'logit'), 
#     data = train)
#   return(model)
# })


# predList <- lapply(modelList, predict, newdata = test, type = "response")
# names(predList) <- str_c("pred_", useTF)

# # predDF <- as_tibble(c(predList, predDSList))
# predDF <- as_tibble(predList)

# test <- bind_cols(test, predDF)



# add model quality and mdoel as tidy DF
byTF <- byTF %>% 
  mutate(tidy_model = map(model, broom::tidy)) %>% 
  mutate(glance = map(model, broom::glance)) %>% 

save(byTF, file = paste0(outPrefix, ".byTF.Rdata"))

# Analyse Model parameters ----------------------------------------------------

# get DF with model quality
modelQualDF <- byTF %>% 
  unnest(glance, .drop = TRUE)
write_tsv(modelQualDF, paste0(outPrefix, ".modelQualDF.tsv"))

# get tidy model DF
modelDF <- byTF %>% 
  unnest(tidy_model, .drop = TRUE)

write_tsv(modelDF, paste0(outPrefix, ".modelDF.tsv"))

# 
# 
# tidyModel <- map(byTF$model, broom::tidy)
# names(tidyModel) <- useTF
# modelDF <- as_tibble(bind_rows(tidyModel, .id = "name"))
# 
# modelDF <- modelDF %>% 
#   mutate(param = str_replace(term, "cor_.*", "cor")) %>% 
#   mutate(param = str_replace(param, "strandOrientation", "")) %>% 
#   mutate(param = parse_factor(param, levels = NULL))
# 
# write_tsv(modelDF, paste0(outPrefix, ".modelDF.tsv"))
modelDF <- modelDF %>% 
  left_join(meta, by = "name")

paramDF <- modelDF %>% 
  group_by(term) %>% 
  summarize(
    n = n(),
    mean = mean(estimate),
    sd = sd(estimate)
  )
write_tsv(paramDF, paste0(outPrefix, ".paramDF.tsv"))

p <- ggplot(modelDF, aes(x = name, y = estimate, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  # geom_text(aes(label = round(estimate, 2)), vjust = "inward") + 
  # facet_grid(param ~ ., scales = "free_y") + 
  geom_text(aes(label = round(estimate, 2)), hjust = "inward") + 
  facet_grid(. ~ term , scales = "free_x") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(outPrefix, ".paramter.barplot.pdf"), w = 6, h = 12)

p <- ggplot(modelDF, aes(x = TF, y = estimate, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = round(estimate, 2)), hjust = "inward") + 
  facet_grid(outType ~ term , scales = "free") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(outPrefix, ".paramter_by_outType.barplot.pdf"), w = 6, h = 12)

# 
# p <- ggplot(modelDF, aes(x = name, y = exp(estimate), fill = name)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   # geom_text(aes(label = round(estimate, 2)), vjust = "inward") + 
#   # facet_grid(param ~ ., scales = "free_y") + 
#   geom_text(aes(label = round(exp(estimate), 2)), hjust = "inward") + 
#   facet_grid(. ~ param , scales = "free_x") + 
#   coord_flip() +
#   labs(y = "Parameter estimate", x = "Model") + 
#   theme_bw() + scale_fill_manual(values = COL_TF) + 
#   theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(p, file = paste0(outPrefix, ".exp_paramter.barplot.pdf"), w = 6, h = 3)


#-------------------------------------------------------------------------------
# Analyse performace for different models -------------------------------------
#-------------------------------------------------------------------------------

# model_names <- c(useTF, paste0(useTF, "_DS"), "dist", "all_TF")
model_names <- useTF

# plot ROC and PRC
modelData <-  mmdata(
  # scores = as.list(select(test, one_of(str_c("pred_", model_names)))),
  scores = map(byTF$data, function(subDF) filter(subDF, fold == 1)[["pred"]]),
  labels = filter(df, fold == 1)[["loop"]],
  modnames = model_names,
  posclass = "Loop"
)

# get rank of models
ranked_models <- auc( evalmod(modelData) ) %>% 
  filter(curvetypes == "PRC") %>% 
  arrange(aucs) %>% 
  select(modnames) %>% 
  unlist()

# build color vector with TFs as names
# TF_models <- ranked_models
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- ranked_models


# build model again with ordered modelnames
modelData <-  mmdata(
  scores = map(byTF$data, function(subDF) filter(subDF, fold == 1)[["pred"]]),
  labels = filter(df, fold == 1)[["loop"]],
  modnames = ranked_models,
  posclass = "Loop"
)

# caluclate ROC, PRC, and basic performance scores
curves <- evalmod(modelData)
# scores <- evalmod(modelData, mode = "basic")

# get AUC of ROC and PRC
aucDF <- as_tibble(auc(curves)) %>%
  left_join(meta, by = c("modnames" = "name")) %>% 
  mutate(modnames = fct_reorder(modnames, aucs, fun = min)) %T>% 
  write_tsv(paste0(outPrefix, ".aucDF.tsv"))

COL_TF <- COL_TF[1:length(model_names)]
names(COL_TF) <- levels(aucDF$modnames)


# %>%
#   mutate(modnames = factor(modnames, ranked_models))

# barplot of AUCs of ROC and PRC

p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)

p <- ggplot(aucDF, aes(x = TF, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ outType, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)


# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curves, "ROC") + 
  # coord_cartesian(expand = FALSE) + 
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFroc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) + 
  theme(legend.position=c(.75,.4))
# g
ggsave(g, file= paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

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

