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
require(colorRamps)  # to pick different colors
require(scales)     # for scales in plotting and percent() formatin function
require(stringr)    # for string function and regular expressions

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
FACTORBOOK_MOTIFS <- TRUE
CTCF_motif_file <- "data/factorbook/factorbookMotifPos.txt.CTCF.bed"
MIN_MOTIF_SCORE <- 3

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
# pie(rep(1, length(COL_TF)), col=COL_TF, labels=COL_TF, main=length(COL_TF))

# COL_TF = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# str_view(grDevices::colors(), "[^\d]$")
# str_view(grDevices::colors(), "^((?!gr(a|e)y).)*[^\\d]$")
# COL_TF <- grDevices::colors()[str_detect(grDevices::colors(), "^((?!(gr(a|e)y|white)).)*[^\\d]$")]
# COL_TF <- COL_TF[!COL_TF %in% c("white", "black", "aliceblue", "azure", "beige", "bisque", "cornsilk", "cyan", "darkorchid", "coral", "darkmagenta")]
# pie(rep(1, 20), col=COL_TF, labels=COL_TF)
# pie(rep(1, length(COL_TF)), col=COL_TF)
# COL_TF = colorRamps::primary.colors(12)

#barplot(1:10, col=COL_TF)
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path(
  "results", 
  paste0("v01_UCSC", 
         ifelse(
           FACTORBOOK_MOTIFS, 
           paste0("_factorbook_", MIN_MOTIF_SCORE), 
           "")))

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

ucscMetaFile <- "data/ENCODE_UCSC_FILES.tsv"


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

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  
  if (!FACTORBOOK_MOTIFS) {
  
    ancGR <- chromloop::motif.hg19.CTCF
    
    # # filter for p-valu <= 10^-6
    # ancGR <- ancGR[ancGR$sig >= 6]
  
    seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
  
  } else {
  
    seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
    ancGR <- import(CTCF_motif_file, format="BED", seqinfo=seqInfo)
    
    # plot motif score distribution    
    p <- tibble(score = score(ancGR)) %>% 
      ggplot(aes(x = score)) + 
      geom_histogram()
    ggsave(p, file = paste0(outPrefix, ".motif_score.histogram.pdf"), w = 7, h = 7)
    
  
    # filter for high scoring motifs
    ancGR <- ancGR[score(ancGR) >= MIN_MOTIF_SCORE]
  
    # flip strand (to match motif of Rao et al. 2014)
    strand(ancGR) <- ifelse(strand(ancGR) == "+", "-", "+")
    ancGR <- sort(ancGR)
  
  }  
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
                              seqinfo = seqInfo)
  )
  
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
  
  for (i in seq_len(nrow(ucscMeta))) {
  
    stopifnot(file.exists(ucscMeta$filePath[i]))
    
    regions(gi) <- chromloop::addCovToGR(
      regions(gi), 
      ucscMeta$filePath[i], 
      window = 1000, 
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
  
  # Annotae with correlation across TFs -------------------------------------
  
  covCols <- paste0("cov_", ucscMeta$TF)
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
    id = row_number(),
    # loopLgt = Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop"
    loop = factor(
      # Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop" | Loop_Mifsud2015_GM12878 == "Loop",
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop"))
  ) %>% 
  # select(id, loop, everything(), -Loop_Rao_GM12878, -Loop_Tang2015_GM12878)
  select(id, loop, everything())

# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = TF, value = cor) 

tidyDF <- tidyDF %>% 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(TF = str_sub(TF, 5))


#-------------------------------------------------------------------------------
# analyse number of positives and negatives -------------------
#-------------------------------------------------------------------------------

vennList <- map(
  select(df, starts_with("Loop_")),
  function(x) which(x == "Loop")
)

venn.diagram(
  x = vennList, 
  filename = paste0(outPrefix, ".loop_balance.venn.png"),
  imagetype = "png",
  euler.d = TRUE, scaled = TRUE
)

vennMat <- df %>% 
  select(starts_with("Loop_")) %>% 
  mutate_all( function(x) x == "Loop")

pdf(paste0(outPrefix, ".loop_balance.venneuler.pdf"))
plot(venneuler(vennMat))
dev.off()

# pie chart of percent positives

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

p <- ggplot(count(df, loop), aes( x ="", y = n, fill = loop)) +
  geom_bar(stat = "identity", width = 0.5) + 
  coord_polar("y") +
  scale_fill_manual(values = COL_LOOP) + 
  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = cumsum(rev(n) / 2),
                label = paste(rev(n), "\n", percent(rev(n) / nrow(df)))),
            size = 5)

ggsave(p, file = paste0(outPrefix, ".loop_balance.pie.pdf"), w = 7, h = 7)


# ------------------------------------------------------------------------------
# NAs in correlation ------------------------------------------------------
# ------------------------------------------------------------------------------

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
# p
ggsave(p, file = paste0(outPrefix, ".percentNotNA.by_TF_and_loop.barplot.pdf"), w = 7, h = 7)

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

p <- ggplot(tidySubDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = TF), lwd = 1.5) + 
  geom_boxplot(fill = "white", lwd = 1.5, width = .2) +
  facet_grid(. ~ TF) + 
  scale_fill_manual(values = COL_TF, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size=20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Correlation of ChIP-seq signal", x="")
  # geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
# p
ggsave(p, file = paste0(outPrefix, ".cor.by_TF_and_loop.boxplot.pdf"), w = 14, h = 7)

# Predict loops --------------------------------------------------------

# useTF <- ucscMeta$name
useTF <- c(ucscMeta$name, "across_TFs")
# useTF <- c("CTCF", "Stat1")


# take same number of true and false interactions for training and validation
nLoop <- sum(df$loop == "Loop")

# downsample negative set to the same size as positves
# predDF <- df %>%
#   group_by(loop) %>%
#   sample_n(nLoop) %>%
#   ungroup()

predDF <- df 
# %>% 
#   filter(strandOrientation == "convergent")



# dived in training and test set by taking 9/10 for training and 1/10 for testing
testCases <- sample.int(nrow(predDF), size = round(nrow(predDF)/10))
train <- predDF[-testCases,]
test <- predDF[testCases,]

# get predictions from single TF with dist

# for (name in ucscMeta$name) {
for (name in useTF) {
  
  message("INFO: Fit model for TF: ", name)
    
  # design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", " dist + cor_", name))
  
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + cor_", name))
  
  model <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = train)

  message("INFO: Predict loops with model: ", name)
  
  test[, str_c("pred_", name)] <- predict(model, newdata = test, type = 'response')
  
    # add_predictions(model, var = str_c("pred_", name))

  # classical prediction
  # fitted.results <- predict(model, newdata=test[,c("dist", str_c("cor_", name))], type = 'response')
  # pred <- ifelse(fitted.results > 0.5,1,0)
  # 
  # test[, str_c("pred_", name)] <- fitted.results
}

# add modles of all TF and only dist

designAll <- as.formula(paste0("loop ~ dist + strandOrientation + ", paste(paste0("cor_", useTF), collapse = " + ")))
modelAll <- glm(
  designAll,
  family = binomial(link = 'logit'),
  data = train)

designDist <- as.formula(paste0("loop ~ dist + strandOrientation"))
# designDist <- as.formula(paste0("loop ~ dist"))
modelDist <- glm(
  designDist, 
  family = binomial(link = 'logit'), 
  data = train)

test[, "pred_all_TF"] <- predict(modelAll, newdata = test, type = 'response')
test[, "pred_dist"] <- predict(modelDist, newdata = test, type = 'response')

# test <- test %>% 
#   mutate(
#     pred_all <- predict(modelAll, newdata = ., type = 'response'),
#     pred_dist <- predict(modelDist, newdata = ., type = 'response')
#   )

#-------------------------------------------------------------------------------
# Analyse performace for different models -------------------------------------
#-------------------------------------------------------------------------------

model_names <- c(useTF, "dist", "all_TF")

# plot ROC and PRC
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

# build color vector with TFs as names
non_TF_models <- c("dist", "all_TF", "across_TFs")
TF_models <- ranked_models[!ranked_models %in% non_TF_models]
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- c(TF_models, non_TF_models)


# build model again with ordered modelnames
modelData <-  mmdata(
  scores = as.list(select(test, one_of(str_c("pred_", ranked_models)))),
  labels = test$loop,
  modnames = ranked_models,
  posclass = "Loop"
)

# caluclate ROC and PRC
curves <- evalmod(modelData)

# get AUC of ROC and PRC
aucDF <- auc(curves) %>%
  left_join(ucscMeta, by = c("modnames" = "name")) %>% 
  mutate(modnames = factor(modnames, ranked_models))

# barplot of AUCs of ROC and PRC

p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = round(aucs, 2)), vjust = 1.5) + 
  facet_grid(curvetypes ~ ., scales = "free_x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_fill_manual(values = COL_TF)
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 7, h = 7)


# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

data(P10N10)
sscurves <- evalmod(scores = P10N10$scores, labels = P10N10$labels)
autoplot(sscurves, "ROC") + 
  theme(legend.position=c(.5,.25))

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

g <- autoplot(curves, "PRC", size=4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFprc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position="none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w=5, h=5)

