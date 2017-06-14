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
require(venneuler)  # for Venn-Euler diagram (area propotional)
require(VennDiagram)# for VennDiagrams
require(BiocParallel) # for parallelisation
require(biobroom)     # to make BioC classes tidy
require(grid)         # for textGrob

#-----------------------------------------------------------------------
# Options for parallel computation
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed=34312)
# set options
register(multicorParam)
# bpparam() # to print current options


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
FACTORBOOK_MOTIFS <- FALSE
CTCF_motif_file <- "data/factorbook/factorbookMotifPos.txt.CTCF.bed"
MIN_MOTIF_SCORE <- 2
MIN_MOTIF_SIG <- 6
OUTPUT_TYPES = FALSE


WINDOW_SIZE <- 1000
BIN_SIZE <- 1
TRUE_LOOPS <- "HIC_ChIAPET"
# TRUE_LOOPS <- "HIC_ChIAPET_CaptureC"

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "green", "darkgreen", "gray30")

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
names(COL_LOOP) <- c("No loop", "Loop")

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path(
  "results", 
  paste0("v02_",
          ifelse(OUTPUT_TYPES, "outTypes", "UCSC"), 
          ifelse(
           FACTORBOOK_MOTIFS, 
           paste0("_factorbook_", MIN_MOTIF_SCORE),
           ifelse(MIN_MOTIF_SIG > 5, paste0("_motifSig", MIN_MOTIF_SIG), "")
          ),
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

ucscMetaFile <- "data/ENCODE_UCSC_FILES.tsv"
outTypeMetaFile <- "data/ENCODE/metadata.fltOuttype.tsv"
# metaFile <- ifelse(OUTPUT_TYPES, outTypeMetaFile, ucscMetaFile)

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

if (OUTPUT_TYPES) {
  # parse metatdata table for input files
  meta <- read_tsv(outTypeMetaFile,
                          col_types = cols(
                            `File accession` = col_character(),
                            TF = col_character(),
                            `Output type` = col_character(),
                            file_nrep = col_character(),
                            exp_nrep = col_integer(),
                            Lab = col_character()
                          )
  )
  
  # add name 
  meta <- meta %>% 
    mutate(name = str_c(TF, "_", str_replace_all(`Output type`, "[ -]", "_"))) %>% 
    select(name, everything())
  
} else {
  # parse ucscMeta file
  meta <- read_tsv(ucscMetaFile,
                   col_types = cols(
                     file = col_character(),
                     TF = col_character()
                   ))
  
  meta <- meta %>% 
    mutate(filePath = file.path("data", "ENCODE", "UCSC", file)) %>% 
    mutate(name = TF)
}

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  
  if (!FACTORBOOK_MOTIFS) {
  
    ancGR <- chromloop::motif.hg19.CTCF

    # plot motif score distribution
    sigDF <- tibble(sig = ancGR$sig)
    p <- sigDF %>% 
      ggplot(aes(x = sig)) + 
      geom_histogram() + 
      geom_vline(xintercept = MIN_MOTIF_SIG) +
      labs(x="Motif hit significance [-log10(p-value)]") + 
      annotation_custom(textGrob("there"), 
                        xmin = MIN_MOTIF_SIG, xmax = Inf, ymin = -Inf, ymax = Inf)
    ggsave(p, file = paste0(outPrefix, ".ancGR.motif_sig.histogram.pdf"), w = 7, h = 7)
    
    sigDF %>% 
      count(sig >= MIN_MOTIF_SIG) %>% 
      write_tsv(path = paste0(outPrefix, ".ancGR.motif_sig.n.tsv"))
        
    # filter for p-valu <= MIN_MOTIF_SIG
    ancGR <- ancGR[ancGR$sig >= MIN_MOTIF_SIG]
  
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

  # add motif score
  gi <- chromloop::addMotifScore(gi, colname = "sig")
  
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  
  trueLoopsTang2015_CTCF <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[1]], seqinfo = seqInfo)
  trueLoopsTang2015_PolII <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[2]], seqinfo = seqInfo)
  trueLoopsTang2015 <- c(trueLoopsTang2015_CTCF, trueLoopsTang2015_PolII)
  
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
  gi <- addInteractionSupport(gi, trueLoopsTang2015_CTCF, "Loop_Tang2015_GM12878_CTCF")
  gi <- addInteractionSupport(gi, trueLoopsTang2015_PolII, "Loop_Tang2015_GM12878_PolII")
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
    
    message("INFO: --> Working on sample: ", meta$name[i], " <--")
  
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
    id = row_number(),
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
  select(id, loop, everything())

# save(df, file = paste0(outPrefix, ".df.Rdata"))  
# load(paste0(outPrefix, ".df.Rdata"))  


# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = TF, value = cor) 

tidyDF <- tidyDF %>% 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(TF = str_sub(TF, 5))


#-------------------------------------------------------------------------------
# analyse number of positives and negatives -------------------
#-------------------------------------------------------------------------------

# simple stats:
countsDF <- tibble(
  nMotif = length(regions(gi)),
  nPairs = nrow(df),
  HiC = sum(df$Loop_Rao_GM12878 == "Loop"),
  HiC_percent = HiC / nPairs * 100, 
  HIC_ChIAPET = sum(df$HIC_ChIAPET == "Loop"),
  HIC_ChIAPET_percent = HIC_ChIAPET / nPairs * 100,
  HIC_ChIAPET_CaptureC = sum(df$HIC_ChIAPET_CaptureC == "Loop"),
  HIC_ChIAPET_CaptureC_percent = HIC_ChIAPET_CaptureC / nPairs * 100
)

write_tsv(countsDF, path = paste0(outPrefix, ".countsDF.tsv"))

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


# Venn-Diagram of ChIA-PET and Hi-C loops --------------------------------------

vennList <- map(
  select(df, Loop_Rao_GM12878, Loop_Tang2015_GM12878_CTCF, Loop_Tang2015_GM12878_PolII),
  function(x) which(x == "Loop")
)
names(vennList) <- c("Hi-C", "CTCF ChIA-PET", "PolII ChIA-PET")

venn.diagram(
  x = vennList, 
  filename = paste0(outPrefix, ".loop_balance.Hi-C_ChIA-PET.venn.png"),
  imagetype = "png",
  euler.d = TRUE, scaled = TRUE,
  lwd = 6,
  lty = 'blank',
  fill = c("cornflowerblue", "green", "yellow"),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1,
  cat.just = list(c(0.6,1) , c(.7,1.2) , c(1,0))
)

vennMat <- df %>% 
  select(Loop_Rao_GM12878, Loop_Tang2015_GM12878_CTCF, Loop_Tang2015_GM12878_PolII) %>% 
  mutate_all( function(x) x == "Loop")
names(vennMat) <- c("Hi-C", "CTCF ChIA-PET", "PolII ChIA-PET")
  
pdf(paste0(outPrefix, ".loop_balance.Hi-C_ChIA-PET.venneuler.pdf"))
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

nDF <- count(df, loop)
p <- ggplot(nDF, aes( x ="", y = n, fill = loop)) +
  geom_bar(stat = "identity", width = 0.5) + 
  coord_polar("y") +
  scale_fill_manual(values = COL_LOOP) + 
  geom_text(aes(y = cumsum(rev(n) / 2),
                label = paste(rev(loop), "\n", rev(n), "\n", percent(rev(n) / nrow(df)))),
            size = 5) + 
  geom_text(aes(x = 0.5, y = 0, label = paste("n =", nrow(df)))) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank()
  )

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
# Analyse motif scores in looping and non looping regions
# ------------------------------------------------------------------------------
# label regions if partisipate in loop
loopAnchros <- anchors(gi[df$loop == "Loop"], id = TRUE, type = "both")
loopAnchrosIds <- unique(unlist(loopAnchros))

regions(gi)$loop_anchor <- seq_along(regions(gi)) %in% loopAnchrosIds
regions(gi)$loop_anchor <- factor(regions(gi)$loop_anchor, c(TRUE, FALSE), c("Loop", "No loop"))

ancDF <- broom::tidy(regions(gi)) 
# %>% 
#   mutate(loop_anchor = factor(loop_anchor, c(TRUE, FALSE), c("Loop", "No loop")))

p <- ggplot(ancDF, aes(x = loop_anchor, y = sig, color = loop_anchor)) + 
  geom_boxplot() + 
  theme_bw() + scale_color_manual(values = COL_LOOP) +
  labs(y = "Significane of Motif hit [-log10(p)]")
ggsave(p, file = paste0(outPrefix, ".motif_signif_by_loop.boxplot.pdf"), w = 3.5, h = 7)

min_anchor_value <- function(gi, colName){
  
  score <- mcols(regions(gi))[,colName]
  a1 <- anchors(gi, id=TRUE, type="first")
  a2 <- anchors(gi, id=TRUE, type="second")
  
  min_score <- pmin(score[a1], score[a2])  
  return(min_score)
}
  
df$min_sig <- min_anchor_value(gi, "sig")

p <- ggplot(df, aes(x = min_sig, color = loop)) + 
  geom_density() + 
  xlim(5, 8) + 
  theme_bw() + scale_color_manual(values = COL_LOOP) +
  labs(x = "Significane of min motif hit [-log10(p)]")
ggsave(p, file = paste0(outPrefix, ".gi.min_motif_th_by_loop.pdf"), w = 6, h = 3)

loop <- df$loop == "Loop"

balanceDF <- tibble(
  min_sig = seq_range(c(5,8), 9, pretty = TRUE),
  n = map_int(min_sig, function(x) sum(df$min_sig >= x)),
  nLoop = map_int(min_sig, function(x) sum(df$min_sig >= x & loop)),
  percent = nLoop / n * 100
)

p <- ggplot(balanceDF, aes(x = min_sig, y = percent), ) + 
  geom_line(color = COL_LOOP[2]) +
  geom_point(color = COL_LOOP[2]) + 
  theme_bw() + labs(
    x = "Motif threshold [-log10(p-value)]",
    y = "Percent positives")
ggsave(p, file = paste0(outPrefix, ".percent_positives_by_moitf_th.pdf"), w = 6, h = 3)

pn <- ggplot(balanceDF, aes(x = min_sig, y = n, label = n)) + 
  geom_line() + 
  geom_point() +
  geom_text(hjust = "left", vjust = "bottom") + 
  theme_bw() + labs(
    x = "Motif threshold [-log10(p-value)]",
    y = "Number of pairs")
ggsave(pn, file = paste0(outPrefix, ".n_by_moitf_th.pdf"), w = 6, h = 3)

pp <- ggplot(balanceDF, aes(x = min_sig, y = nLoop, label = nLoop)) + 
  geom_line(color = COL_LOOP[2]) + 
  geom_point(color = COL_LOOP[2]) +
  geom_text(hjust = "left", vjust = "bottom", color = COL_LOOP[2]) + 
  theme_bw() + labs(
    x = "Motif threshold [-log10(p-value)]",
    y = "Number of positive pairs")
ggsave(pp, file = paste0(outPrefix, ".positives_by_moitf_th.pdf"), w = 6, h = 3)

# plot_grid(p, pn, labels=c("A", "B"), ncol = 1, nrow = 2)
# grid.arrange(p, pn, ncol = 1, nrow = 2)

# strand orientation 

p <- df %>%
  group_by(loop, strandOrientation) %>%
  summarise(
    n = n()
  ) %>% 
  ggplot(aes(x = strandOrientation, y = n, fill = strandOrientation)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = n), vjust = 2) +
  facet_grid(loop ~ . , scales = "free_y") + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") + 
  labs(x = "Motif orientation", y = "Number of pairs")

ggsave(p, file = paste0(outPrefix, "n_by_strand_and_loop.barplot.pdf"),
       w = 5, h = 5)

# distance 
p <- ggplot(df, aes(x = dist / 1000, fill = loop, color = loop)) +
  geom_density(size = 2, alpha = 0.5) +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = "bottom") + 
  labs(y = "Density", x = "Loop Distance [kb]") 
ggsave(p, file = paste0(outPrefix, ".distance_by_loop.density.pdf"), w = 6, h = 3)

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

# ------------------------------------------------------------------------------
# Predict loops --------------------------------------------------------
# ------------------------------------------------------------------------------

useTF <- c(meta$name, "across_TFs")
# useTF <- c("CTCF", "Stat1")


# dived in training and test set by taking 9/10 for training and 1/10 for testing
testCases <- sample.int(nrow(df), size = round(nrow(df)/10))
train <- df[-testCases,]
test <- df[testCases,]

# nLoop <- sum(train$loop == "Loop")
# # downsample negative set to the same size as positves
# trainDS <- train %>% 
#     group_by(loop) %>%
#     sample_n(nLoop) %>%
#     ungroup()
  
# get predictions from single TF with dist

# for (name in ucscMeta$name) {
# for (name in useTF) {
# predList <- bplapply(useTF, function(name) {
predList <- bplapply(useTF, function(name) {
    
  message("INFO: Fit model for TF: ", name)
    
  # design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", " dist + cor_", name))
  
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + score_min + cor_", name))
  
  model <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = train)
  
  message("INFO: Predict loops with model: ", name)
  
  # test[, str_c("pred_", name)] <- predict(model, newdata = test, type = 'response')
  return( predict(model, newdata = test, type = 'response') )
})
names(predList) <- str_c("pred_", useTF)

# predDF <- as_tibble(c(predList, predDSList))
predDF <- as_tibble(predList)

test <- bind_cols(test, predDF)

# add modles of all TF and only dist

designAll <- as.formula(paste0("loop ~ dist + strandOrientation + ", paste(paste0("cor_", useTF), collapse = " + ")))
modelAll <- glm(
  designAll,
  family = binomial(link = 'logit'),
  data = train)
test[, "pred_all_TF"] <- predict(modelAll, newdata = test, type = 'response')

modelDistOrientation <- glm(
  loop ~ dist + strandOrientation + score_min, 
  family = binomial(link = 'logit'), 
  data = train)
test[, "pred_Dist+Orientation+Motif"] <- predict(modelDistOrientation, newdata = test, type = 'response')

modelDistOnly <- glm(
  loop ~ dist, 
  family = binomial(link = 'logit'), 
  data = train)
test[, "pred_Dist"] <- predict(modelDistOnly, newdata = test, type = 'response')

modelOrientation <- glm(
  loop ~ strandOrientation, 
  family = binomial(link = 'logit'), 
  data = train)
test[, "pred_Orientation"] <- predict(modelOrientation, newdata = test, type = 'response')

modelMotif <- glm(
  loop ~ score_min, 
  family = binomial(link = 'logit'), 
  data = train)
test[, "pred_Motif"] <- predict(modelOrientation, newdata = test, type = 'response')

save(train, test, file = paste0(outPrefix, ".train_test.Rdata"))

#-------------------------------------------------------------------------------
# Model parameters
#-------------------------------------------------------------------------------
modelList <- bplapply(useTF, function(name) {
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + score_min + cor_", name))
  
  model <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = train)
  return(model)
})

tidyModel <- lapply(modelList, broom::tidy)
names(tidyModel) <- useTF
modelDF <- as_tibble(bind_rows(tidyModel, .id = "name"))

modelDF <- modelDF %>% 
  mutate(param = str_replace(term, "cor_.*", "cor")) %>% 
  mutate(param = str_replace(param, "strandOrientation", "")) %>% 
  mutate(param = parse_factor(param, levels = NULL))

p <- ggplot(modelDF, aes(x = name, y = estimate, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  # geom_text(aes(label = round(estimate, 2)), vjust = "inward") + 
  # facet_grid(param ~ ., scales = "free_y") + 
  geom_text(aes(label = round(estimate, 2)), hjust = "inward") + 
  facet_grid(. ~ param , scales = "free_x") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(outPrefix, ".paramter.barplot.pdf"), w = 6, h = 3)
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
model_names <- c(useTF, "all_TF", "Dist+Orientation", "Dist", "Orientation")

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

COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray30", "lightblue", "orange", "lightgreen", "lightred", "gray70")

# build color vector with TFs as names
non_TF_models <- c("all_TF", "Dist+Orientation+Motif", "Dist", "Orientation", "Motif", "across_TFs")
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

# caluclate ROC, PRC, and basic performance scores
curves <- evalmod(modelData)
# scores <- evalmod(modelData, mode = "basic")

# get AUC of ROC and PRC
aucDF <- auc(curves) %>%
  left_join(meta, by = c("modnames" = "name")) %>%
  mutate(modnames = factor(modnames, ranked_models))

# barplot of AUCs of ROC and PRC

p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), vjust = 1.5) +
  facet_grid(curvetypes ~ ., scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot_new.pdf"), w = 7, h = 7)


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
                     guide = guide_legend(override.aes = list(size = 1.5),
                                          reverse = TRUE)) + 
  theme(legend.position=c(.75,.45)) # , legend.text = element_text(size = 7)
# g
ggsave(g, file= paste0(outPrefix, ".ROC_new.pdf"), w = 5, h = 5)

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
ggsave(g, file = paste0(outPrefix, ".PRC_new.pdf"), w = 5, h = 5)

