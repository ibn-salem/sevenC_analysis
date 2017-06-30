################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(rtracklayer)  # to import() BED files
require(BiocParallel) # for parallelisation

require(RColorBrewer)   # for nice colors
require(colorRamps)  # to pick different colors
require(precrec)      # for ROC and PRC curves
require(scales)     # for scales in plotting and percent() formatin function

# require(venneuler)  # for Venn-Euler diagram (area propotional)
# require(VennDiagram)# for VennDiagrams

# require(biobroom)     # to make BioC classes tidy
require(grid)         # for textGrob
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling

#-------------------------------------------------------------------------------
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
TRUE_LOOPS <- "HIC_ChIAPET"

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path("results", 
  paste0(
    "input_types_v01",
    ifelse(MIN_MOTIF_SIG > 5, paste0("_motifSig", MIN_MOTIF_SIG), ""),
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


DATA_TYPES_META_FILE = "data/DATA_TYPES_metadata.tsv"

#-------------------------------------------------------------------------------
# Parse and filter meta data  -----------------------------------
#-------------------------------------------------------------------------------

meta <- read_tsv(DATA_TYPES_META_FILE, col_types = cols(
  file_accession = col_character(),
  data_type = col_character(),
  TF = col_character(),
  output_type = col_character(),
  path = col_character()
))
                  
meta <- meta %>% 
  mutate(
    name = str_replace_all(str_c(data_type, TF, output_type, sep = "_"), "[ -]", "_"),
    file_exists = file.exists(path)
    ) %>% 
  select(name, file_exists, everything()) %>% 
  filter(file_exists)

write_tsv(meta, paste(outPrefix, ".meta_filtered.tsv"))

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  ancGR <- chromloop::motif.hg19.CTCF

  # filter for p-valu <= MIN_MOTIF_SIG
  ancGR <- ancGR[ancGR$sig >= MIN_MOTIF_SIG]

  seqInfo <- seqinfo(ancGR)
  

  # get all pairs within 1M distance
  gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)
  
  # add strand combinations
  gi <- chromloop::addStrandCombination(gi)

  # add motif score
  gi <- chromloop::addMotifScore(gi, colname = "sig")
  
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  
  trueLoopsTang2015_CTCF <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[1]], seqinfo = seqInfo)
  trueLoopsTang2015_PolII <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[2]], seqinfo = seqInfo)
  trueLoopsTang2015 <- c(trueLoopsTang2015_CTCF, trueLoopsTang2015_PolII)
  
  gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
  gi <- addInteractionSupport(gi, trueLoopsTang2015_CTCF, "Loop_Tang2015_GM12878_CTCF")
  gi <- addInteractionSupport(gi, trueLoopsTang2015_PolII, "Loop_Tang2015_GM12878_PolII")
  gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))
  
}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Annotae with coverage and correlation -------------------------------------

if (!GI_LOCAL ) {
  
  for (i in seq_len(nrow(meta))) {
    
    message("INFO: --> Working on sample: ", meta$name[i], " <--")

    stopifnot(file.exists(meta$path[i]))
    
    regions(gi) <- chromloop::addCovToGR(
      regions(gi), 
      meta$path[i], 
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
    id = row_number(),
    loop = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop"))
  ) %>% 
  select(id, loop, everything())

save(df, file = paste0(outPrefix, ".df.Rdata"))  
# load(paste0(outPrefix, ".df.Rdata"))

# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = name, value = cor) 

tidyDF <- tidyDF %>% 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(name = str_sub(name, 5))


# ------------------------------------------------------------------------------
# Compare Correlation with Boxplot
# ------------------------------------------------------------------------------

tidySubDF <- tidyDF %>% 
  sample_n(min(nrow(tidyDF), 10^7))

p <- ggplot(tidySubDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = name), lwd = 1.5) + 
  geom_boxplot(fill = "white", lwd = 1.5, width = .2) +
  facet_grid(. ~ name) + 
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(12, "Set3"))(nrow(meta)), 
    guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size=20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Correlation of ChIP-seq signal", x="")
  # geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
# p
ggsave(p, file = paste0(outPrefix, ".cor.by_name_and_loop.boxplot.pdf"), w = 14, h = 7)

# ------------------------------------------------------------------------------
# Predict loops --------------------------------------------------------
# ------------------------------------------------------------------------------


# dived in training and test set by taking 9/10 for training and 1/10 for testing
testCases <- sample.int(nrow(df), size = round(nrow(df)/10))
train <- df[-testCases,]
test <- df[testCases,]


predList <- lapply(meta$name, function(name) {
    
  message("INFO: Fit model for TF: ", name)
    
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + score_min + cor_", name))
  
  model <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = train)
  
  message("INFO: Predict loops with model: ", name)
  
  return( predict(model, newdata = test, type = 'response') )
})

names(predList) <- str_c("pred_", meta$name)

test <- bind_cols(test, as_tibble(predList))

# add modles of all TF and only dist
designList <- list(
  "Dist+Orientation+Motif" = loop ~ dist + strandOrientation + score_min,
  "Dist" =  loop ~ dist,
  "Orientation" = loop ~ strandOrientation,
  "Motif" = loop ~ score_min
)

for (modName in names(designList)) {
  design <- designList[[modName]]
  model <- glm(
    design,
    family = binomial(link = 'logit'),
    data = train)
  test[, str_c("pred_", modName)] <- predict(model, newdata = test, type = 'response')
}

################################################################################
save(train, test, file = paste0(outPrefix, ".train_test.Rdata"))
# load(paste0(outPrefix, ".train_test.Rdata"))
# useTF <- meta$name
################################################################################

#-------------------------------------------------------------------------------
# Analyse performace for different models -------------------------------------
#-------------------------------------------------------------------------------

model_names <- c(meta$name, names(designList))

# get mmdata object
modelData <-  mmdata(
  scores = as.list(select(test, one_of(str_c("pred_", model_names)))),
  labels = test$loop,
  modnames = model_names,
  posclass = "Loop"
)

# get rank of models
ranked_models <- auc( evalmod(modelData) ) %>% 
  filter(curvetypes == "PRC") %>% 
  arrange(aucs) %>% 
  select(modnames) %>% 
  unlist()

# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray30", "cornflowerblue", "orange", "lightgreen", "indianred1", "gray70")
# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray30", "cornflowerblue", "orange", "lightgreen", "indianred1")
COL_TF <- colorRampPalette(brewer.pal(12, "Set3"))(length(model_names))

# build color vector with TFs as names
# non_TF_models <- c("all_TF", "Dist+Orientation+Motif", "Dist", "Orientation", "Motif", "across_TFs")
# non_TF_models <- c("all_TF", "Dist+Orientation+Motif", "Dist", "Orientation", "Motif")
TF_models <- ranked_models[!ranked_models %in% names(designList)]
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- c(TF_models, names(designList))


# build model again with ordered modelnames
modelData <-  mmdata(
  scores = as.list(dplyr::select(test, one_of(str_c("pred_", ranked_models)))),
  labels = test$loop,
  modnames = ranked_models,
  posclass = "Loop"
)

# caluclate ROC, PRC, and basic performance scores
curves <- evalmod(modelData)
# scores <- evalmod(modelData, mode = "basic")

# get AUC of ROC and PRC
aucDF <- as_tibble(auc(curves)) %>%
  left_join(meta, by = c("modnames" = "name")) %>%
  mutate(modnames = factor(modnames, ranked_models))

# barplot of AUCs of ROC and PRC
p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), vjust = 1.5) +
  facet_grid(curvetypes ~ ., scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot_new.pdf"), w = 7, h = 7)

# barplot of AUCs by output type
p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = TF)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), vjust = 1.5) +
  facet_grid(curvetypes ~ output_type, scales = "free", space = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        legend.position = "right") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Models", y = "AUC")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot_by_data_type.pdf"), w = 14, h = 7)


# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curves, "ROC") + 
  # coord_cartesian(expand = FALSE) + 
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       " (AUC=", signif(aucDFroc$aucs,3), ")"),
                     guide = guide_legend(override.aes = list(size = 1.5),
                                          reverse = TRUE)) + 
  theme(legend.position = c(.75,.45)) # , legend.text = element_text(size = 7)
#g
ggsave(g, file = paste0(outPrefix, ".ROC_new.pdf"), w = 5, h = 5)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(curves, "PRC", size = 4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFprc$modnames, 
                       " (AUC=", signif(aucDFprc$aucs,3)), ")",
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position = "none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC_new.pdf"), w = 5, h = 5)

