################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- TRUE

outPrefix <- file.path("results", "v01")
dir.create(dirname(outPrefix), showWarnings = FALSE)

# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

metadataFile <- "data/ENCODE/metadata.flt.tsv"
# work only on subset here:

useTFs <- c(
  "ZNF143",
  "STAT1",
  "ZNF274",
  "ZNF384",
  "CTCF",
  "RAD21",
  "POLR2A",
  "NFYB"
)

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  ancGR <- chromloop::motif.hg19.CTCF
  
  seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
  
  # get all pairs within 1M distance
  gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)
  
  # add strand combinations
  gi <- chromloop::addStrandCombination(gi)
  
  
  # parse loops
  trueLoops <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  
  # ol <- IRanges::overlapsAny(gi, trueLoops)
  # gi$Loop_Rao_GM12878 <- factor(ol, c(FALSE, TRUE), c("No loop", "Loop"))
  
  gi <- addInteractionSupport(gi, trueLoops, "Loop_Rao_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))

}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Parse and filter input ChiP-seq data  -----------------------------------

# parse metatdata table for input files
meta <- read_tsv(metadataFile,
                 col_types = cols(
                   .default = col_character(),
                   file_nrep = col_integer(),
                   exp_nrep = col_integer(),
                   size = col_double(),
                   `Experiment date released` = col_date(format = ""),
                   `Technical replicate.x` = col_integer(),
                   Size = col_double(),
                   `Biological replicate` = col_number(),
                   `Technical replicate.y` = col_number(),
                   usedTF = col_logical()
                 ))

# # filter input data
# meta <- meta %>% 
#   filter(TF %in% useTFs) %>%
#   mutate(TF = parse_factor(TF, levels = useTFs)) %>%
#   arrange(TF)

# filter input data
meta <- meta %>% 
  filter(TF %in% useTFs) %>%
  filter(`Output type` == "fold change over control") %>%
  # filter(`Output type` == "raw signal") %>%
  mutate(TF = parse_factor(TF, intersect(useTFs, .$TF))) %>%
  arrange(TF)

# # DEBUG use only two examples here:
# meta <- meta %>%
#   filter(TF %in% c("POLR2A", "CTCF"))

# Annotae with coverage and correlation -------------------------------------

for (i in seq_along(meta$TF)) {
  regions(gi) <- chromloop::addCovToGR(
    regions(gi), 
    meta$filePath[i], 
    window = 10, 
    colname = paste0("cov_", meta$TF[i])
    )
}

# compute pairwise correlaitons
for (i in seq_along(meta$TF)) {
  gi <- chromloop::applyToCloseGI(
    gi, 
    datcol = paste0("cov_", meta$TF[i]),
    fun = cor, 
    colname = paste0("cor_", meta$TF[i])
  )
}

# save file for faster reload
if (!GI_LOCAL ) {
  save(gi, file = paste0(outPrefix, ".gi.Rdata"))
} else {
  # load(paste0(outPrefix, ".gi.Rdata"))  
}


# Analyse loopps --------------------------------------------------------

df <- as_tibble(as.data.frame(mcols(gi)))

# plot correlation of CTCF by strand and loop
df %>%
  # sample_n(size = 10^6) %>%
  # filter(Loop_Rao_GM12878 == "Loop") %>%
  ggplot(aes(x = Loop_Rao_GM12878, cor_CTCF, fill = strandOrientation)) + 
  geom_violin() + 
  geom_boxplot(width = 0.2, fill = "white") + 
  facet_grid( ~ strandOrientation, scales = "free_y") + 
  theme_bw()

ggsave(paste0(outPrefix, "cor_CTCF_by_strand_and_Loop_Rao_GM12878.violin.pdf"),
       w = 6, h = 6)

df %>%
  sample_n(size = 10^6) %>%
  # filter(Loop_Rao_GM12878 == "Loop") %>%
  ggplot(aes(dist, fill = strandOrientation)) + 
  geom_histogram(aes(y = ..density..)) + 
  facet_grid(strandOrientation ~ Loop_Rao_GM12878, scales = "free_y") + 
  theme_bw()
ggsave(paste0(outPrefix, "dist_by_strand_and_Loop_Rao_GM12878.histogram.pdf"),
       w = 6, h = 6)

df %>%
  group_by(Loop_Rao_GM12878, strandOrientation) %>%
  summarise(
    n = n()
  ) %>% 
  ggplot(aes(x = strandOrientation, y = n, fill = strandOrientation)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = 2) +
  facet_grid(Loop_Rao_GM12878 ~ . , scales = "free_y") + 
  theme_bw()

ggsave(paste0(outPrefix, "n_by_strand_and_Loop_Rao_GM12878.barplot.pdf"),
       w = 6, h = 6)

# Predict loops --------------------------------------------------------

# take same number of true and false interactions for training and validation
nLoop <- sum(df$Loop_Rao_GM12878 == "Loop")
predDF <- df %>%
  mutate(loop = Loop_Rao_GM12878) %>%
  group_by(loop) %>%
  sample_n(nLoop) 

# dived in training and test set by taking 10.000 randomly for testing
testCases <- sample.int(nrow(predDF), size = round(nrow(predDF)/10))
train <- predDF[-testCases,]
test <- predDF[testCases,]

# # get predictions from single TF with dist
# factPred <- lapply(ORDERED_FACTORS, function(FACT){
FACT <- "CTCF"  

design <- as.formula(paste0("loop ~ ", " dist +", "cor_", FACT))

model <- glm(
  design, 
  family = binomial(link = 'logit'), 
  data = train)

test <- test %>% 
  add_predictions(model)

# Analyse performace --------------------------------------------------

# plot ROC and PRC
modelData <-  mmdata(
  scores = test$pred,
  labels = test$loop,
  modnames = "CTCF + dist",
  posclass = "Loop"
)

# caluclate ROC and PRC
cuves <- evalmod(modelData)

# get AUC of ROC and PRC
aucDF <- auc(cuves)


# plat AUCs
g <- ggplot(aucDF, aes(x = modnames, y = aucs)) +
  geom_bar(stat="identity", color="black") + facet_grid(.~curvetypes) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  geom_text(aes(label=signif(aucs,3)), vjust=1.5)

# ggsave(g, file=paste0(outPrefix, ".", LOOP_COL, ".", CELL, ".ChIP-seq_cov_across_all.cor.ROC_PRC_AUC.barplot.pdf"), w=7, h=7)

# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(cuves, "ROC") + 
  annotate("text", x = .6, y = seq(0, .25, length.out = nrow(aucDFroc)), 
           label=paste0(aucDFroc$modnames, ": AUC=", signif(aucDFroc$aucs,3)))
g

# ggsave(g, file=paste0(outPrefix, ".", LOOP_COL, ".", CELL, ".ChIP-seq_cov_across_all.cor.ROC.curve.pdf"), w=7, h=7)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(cuves, "PRC") + 
  annotate("text", x = .6, y = seq(0, .25, length.out = nrow(aucDFprc)), 
           label = paste0(aucDFprc$modnames, ": AUC=", signif(aucDFprc$aucs,3)))

# ggsave(g, file=paste0(outPrefix, ".", LOOP_COL, ".", CELL, ".ChIP-seq_cov_across_all.cor.PRC.curve.pdf"), w=7, h=7)

