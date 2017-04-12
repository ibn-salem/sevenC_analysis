################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- TRUE

COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
#barplot(1:10, col=COL_TF)
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]

outPrefix <- file.path("results", "v01_UCSC")
dir.create(dirname(outPrefix), showWarnings = FALSE)

# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_GM12878_File <- 
  "data/Tang2015/GSM1872886_GM12878_CTCF_PET_clusters.txt"

ucscMetaFile <- "data/ENCODE_UCSC_FILES.tsv"

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  ancGR <- chromloop::motif.hg19.CTCF
  
  seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
  
  # get all pairs within 1M distance
  gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)
  
  # add strand combinations
  gi <- chromloop::addStrandCombination(gi)
  
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  trueLoopsTang2015 <- chromloop::parseLoopsTang2015(
    LoopTang2015_GM12878_File, seqinfo = seqInfo)
  
  # ol <- IRanges::overlapsAny(gi, trueLoops)
  # gi$Loop_Rao_GM12878 <- factor(ol, c(FALSE, TRUE), c("No loop", "Loop"))
  
  gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
  gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))
  
}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Parse and filter input ChiP-seq data  -----------------------------------

# parse ucscMeta file
ucscMeta <- read_tsv(ucscMetaFile,
                     col_types = cols(
                       file = col_character(),
                       TF = col_character()
                     ))

ucscMeta <- ucscMeta %>% 
  mutate(filePath = file.path("data", "ENCODE", "UCSC", file)) %>% 
  mutate(name = TF)


# Annotae with coverage and correlation -------------------------------------


# add coverage and correalation for UCSC files

# for (i in seq_len(nrow(ucscMeta))) {
for (i in 8:10) {
    
  stopifnot(file.exists(ucscMeta$filePath[i]))
  
  regions(gi) <- chromloop::addCovToGR(
    regions(gi), 
    ucscMeta$filePath[i], 
    window = 10, 
    colname = paste0("cov_", ucscMeta$TF[i])
  )
  
  # add correlations
  gi <- chromloop::applyToCloseGI(
    gi, 
    datcol = paste0("cov_", ucscMeta$TF[i]),
    fun = cor, 
    colname = paste0("cor_", ucscMeta$TF[i])
  )  
  
}

# save file for faster reload
if (!GI_LOCAL ) {
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
    loop = Loop_Tang2015_GM12878
  ) %>% 
  select(id, loop, everything(), -Loop_Rao_GM12878, -Loop_Tang2015_GM12878)


# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = TF, value = cor) 

tidyDF <- tidyDF %>% 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(TF = str_sub(TF, 5))

# NAs in correlation ------------------------------------------------------

naDF <- tidyDF %>% 
  group_by(TF, loop) %>% 
  summarise(
    n = n(),
    nNotNA = sum(!is.na(cor)),
    percentNotNA = nNotNA / n * 100
  ) 

p <- ggplot(naDF, aes(x = TF, y = percentNotNA, fill = TF)) +
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = round(percentNotNA, 2)), vjust = 1.5) +
  facet_grid(loop ~ ., scales = "free_x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_fill_manual(values = COL_TF)
p
ggsave(p, file=paste0(outPrefix, ".percentNotNA.by_TF_and_loop.barplot.pdf"), w = 7, h = 7)


# Compare Correlation with Boxplot
# pvalDF <- ddply(plotDF, .(factor), summarize, p=wilcox.test( as.formula(paste("value ~ ", LOOP_COL)))$p.value)

p <- ggplot(tidyDF, aes(x = loop, y = cor)) +
  geom_boxplot(aes(color=loop), lwd=1.5) + 
  facet_grid(. ~ TF) + 
  scale_color_manual(values=COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size=20), 
    axis.text.x=element_blank(), 
    legend.position = "bottom") + 
  labs(y="Correlation of ChIP-seq signal", x="")
  # geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
p
# ggsave(p, file=paste0(outPrefix, ".", LOOP_COL, ".", CELL, "_w", WINDOW, "_b", BIN_SIZE, ".ChIP-seq_profile.cor.boxplot.pdf"), w=14, h=7)


# Predict loops --------------------------------------------------------

# take same number of true and false interactions for training and validation
nLoop <- sum(df$loop == "Loop")

# downsample negative set to the same size as positves
predDF <- df %>%
  group_by(loop) %>%
  sample_n(nLoop) %>% 
  ungroup() 

# dived in training and test set by taking 9/10 for training and 1/10 for testing
testCases <- sample.int(nrow(predDF), size = round(nrow(predDF)/10))
train <- predDF[-testCases,]
test <- predDF[testCases,]

# get predictions from single TF with dist

for (name in ucscMeta$name) {
  
  # design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + ", "cor_", str_replace(name, "-", "_")))
  design <- as.formula(paste0("loop ~ ", "cor_", str_replace(name, "-", "_")))
  
  model <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = train)
  
  test <- test %>% 
    add_predictions(model, var = str_c("pred_", name))
  
}


# Analyse performace for different output types -------------------------------------

# plot ROC and PRC
modelData <-  mmdata(
  scores = as.list(select(test, one_of(str_c("pred_", ucscMeta$name)))),
  labels = test$loop,
  modnames = ucscMeta$name,
  posclass = "Loop"
)

# caluclate ROC and PRC
cuves <- evalmod(modelData)

# get AUC of ROC and PRC
aucDF <- auc(cuves) %>% 
  left_join(ucscMeta, by = c("modnames" = "name"))

# barplot of AUCs of ROC and PRC

p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = TF)) +
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = round(aucs, 2)), vjust = 1.5) + 
  facet_grid(curvetypes ~ ., scales = "free_x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_fill_manual(values = COL_TF)

ggsave(p, file=paste0(outPrefix, ".AUC_ROC_PRC.by_OutputType_and_TF.barplot.pdf"), w = 7, h = 7)


# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(cuves, "ROC") + 
  annotate("text", x = .6, y = seq(0, .25, length.out = nrow(aucDFroc)), 
           label=paste0(aucDFroc$modnames, ": AUC=", signif(aucDFroc$aucs,3)))

# ggsave(g, file=paste0(outPrefix, ".", LOOP_COL, ".", CELL, ".ChIP-seq_cov_across_all.cor.ROC.curve.pdf"), w=7, h=7)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(cuves, "PRC") + 
  annotate("text", x = .6, y = seq(0, .25, length.out = nrow(aucDFprc)), 
           label = paste0(aucDFprc$modnames, ": AUC=", signif(aucDFprc$aucs,3)))
g


