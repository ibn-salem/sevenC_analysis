################################################################################
# Analysis of predictd chromatin looping interactions by compairing with 
# enhancer-promoter associations
################################################################################


library(chromloop)    # devtools::install_github("ibn-salem/chromloop")
library(rtracklayer)  # to import() BED files
# library(EnsDb.Hsapiens.v75)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for human genes
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(RColorBrewer)   # for nice colors
library(feather)      # for efficient storing of data.frames
library(readxl)       # to read excel files
source("R/chromloop.functions.R")
source("R/gr_associations.R")


# 0) Set parameter --------------------------------------------------------

outPrefix <- "results/v05.TF_enhancer_associations"

# topN_models_file = "results/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.bestNModelDF.tsv"
models_prefix = "results/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1"

selected_model_prefix <- "results/v05_selected_models.motifPval2.5e-06_w1000_b1"
# model_prefix <- "results/v04_screen_TF_lfc.motifSig6_w1000_b1" # bestNModelDF.tsv

# enahncer_association_file <- "data/EnhancerAtlas/GM12878_EP.txt"

enhancer_associations_urls <- "data/EnhancerAtlas_EP_cell-lines_urls.txt"

metadata_file <- "data/ENCODE/metadata.fcDF.tsv"

# --------------- Read predicted loops from selected model analysis ------------
# lood predicted loops wiht coverage annotation
gi <- read_rds(paste0(selected_model_prefix, ".gi.rds"))

# use default model of average across best 10 TFs
model <- read_tsv(paste0(models_prefix, ".bestNModelDF.tsv"))

BestN_cutof <- read_tsv(paste0(models_prefix, ".topNf1ModelDF.tsv")) %>% 
  pull(mean_max_cutoff)

# predict loops using default model
loops <- chromloop::predLoops(
  gi,
  formula = ~ dist + strandOrientation + score_min + cor_RAD21,
  betas = model$estimate_mean,
  cutoff = BestN_cutof
  )

# get regions enclosed by loops
loopRange <- interactionRange(loops)

# extend anchros and find direct interactions
# extended_size = 5*10^4
inner_extended_size = 5000
outer_extended_size = 500
extLoops <- extendAnchors(loops, 
                          inner = inner_extended_size, 
                          outer = outer_extended_size)

# ----------------------- Parse enhancer-gene associations ---------------------
seqInfo <- seqinfo(motif.hg19.CTCF)

epDF <- read_tsv(
    enhancer_associations_urls, 
    col_names = "url", col_types = cols(url = col_character())
  ) %>% 
  mutate(
    file = basename(url), 
    cell = str_replace(file, "_EP.txt", ""),
    path = file.path("data", "EnhancerAtlas", file)
  )

#' A function to read enhancer-gene associations as GInteraction object
readEnhancerAssociations <- function(inFile, seqInfo) {

  ep_col_names <-  c("chromEnh", "chromStart", "chromEnd", "GeneID", 
                     "chromGene", "TSS", "Transcript", "NULL", "signalValue", 
                     "EPscore", "additional_score_1", "additional_score_2", 
                     "additional_score_3", "additional_score_4")  
  epDF <- read_tsv(
    inFile, 
    col_names = ep_col_names,
    col_types = cols(
      chromEnh = col_character(),
      chromStart = col_integer(),
      chromEnd = col_integer(),
      GeneID = col_integer(),
      chromGene = col_character(),
      TSS = col_integer(),
      Transcript = col_character(),
      `NULL` = col_character(),
      signalValue = col_double(),
      EPscore = col_double(),
      additional_score_1 = col_character(),
      additional_score_2 = col_character(),
      additional_score_3 = col_character(),
      additional_score_4 = col_character()
    )
  )
  
  enhancerGR <- GRanges(epDF$chromEnh, IRanges(epDF$chromStart + 1, epDF$chromEnd),
                        seqinfo = seqInfo, score = epDF$signalValue)
  
  tssGR <- GRanges(epDF$chromGene, IRanges(epDF$TSS + 1, epDF$TSS + 1),
                   seqinfo = seqInfo, name = epDF$Transcript, 
                   geneID = epDF$GeneID)
  
  ep <- GInteractions(enhancerGR, tssGR)
  mcols(ep) <- select(epDF, Transcript, signalValue, EPscore)
  
  return(ep)
}

epDF <- epDF %>% 
  mutate(
    ep = map(
      path, 
      readEnhancerAssociations, 
      seqInfo = seqinfo(motif.hg19.CTCF)),
    n = map_int(ep, length)
  ) %>% 
  select(-url, -path)

# filter for at least 30k E-P pairs and top 30k
epDF <- epDF %>% 
  filter(n > 30000) %>% 
  mutate(
    ep_order = map(ep, ~ order(.x$EPscore, decreasing = TRUE)),
    top_idx = map(ep_order, ~ .x[seq(30000)]),
    ep = map2(ep, top_idx, ~ .x[.y]),
    n_flt = map_int(ep, length)
  )

# get range of interactions
epDF <- epDF %>% 
  mutate(
    range = map(ep, interactionRange),
    range_size = map(range, width)
    )

# epRangeList <- map(epList, interactionRange)

#=== Analyse size distribution of enhancers ====================================

p <- ggplot(epDF, aes(x = cell, y = n, fill = cell)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outPrefix, ".E-P_number.barplot.pdf"), w = 6, h = 3)

sizeDF <- epDF %>% 
  select(cell, range_size) %>% 
  unnest(range_size)

p <- ggplot(sizeDF, aes(x = cell, y = range_size, color = cell)) + 
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = "E-P range size (bp)")
ggsave(paste0(outPrefix, ".E-P_range_size.boxplot.pdf"), w = 6, h = 3)

#=== Analyse E-P association with predicted loops ==============================

# hits <- map(epList, findOverlaps, subject = loops)

epDF <- epDF %>% 
  mutate(
    # check if E-P associations are within predicted loops
    rangeHits = map(range, ~findOverlaps(query = .x, subject = loopRange, type = "within")),
    inLoop = map2(rangeHits, ep,  ~ seq(length(.y)) %in% unique(queryHits(.x))), 
    inLoopPercent = map_dbl(inLoop, mean) * 100,
    # analyse percent of Loops that do connect E-Ps
    loopWithEP = map(rangeHits, ~ seq(length(loops)) %in% unique(subjectHits(.x)) ),
    loopWithEPpercent = map_dbl(loopWithEP, mean) * 100
    )


epDF <- epDF %>% 
  mutate(
    # check if E-P associations overlap directly with anchor extended loops
    pairdOvlHits = map(ep, findOverlaps, subject = extLoops, ignore.strand = FALSE),
    inExtLoop = map2(pairdOvlHits, ep,  ~ seq(length(.y)) %in% unique(queryHits(.x))), 
    inExtLoopPercent = map_dbl(inExtLoop, mean) * 100,
    # analyse percent of extended loops that do connect E-Ps
    extLoopWithEP = map(pairdOvlHits, ~ seq(length(extLoops)) %in% unique(subjectHits(.x)) ),
    extLoopWithEPpercent = map_dbl(extLoopWithEP, mean) * 100
  )

# rangeHits <- map(epRangeList, findOverlaps, subject = loopRange, type = "within")

# inLoop <- map2(rangeHits, epList, ~ seq(length(.y)) %in% unique(subjectHits(.x)) )

# inLoopPercent <- map(inLoop, mean)

# loopWithEP <- map(rangeHits, ~ seq(length(loops)) %in% unique(queryHits(.x)) )
# loopWithEPpercent <- map(loopWithEP, mean)

# -------------- plot percent of E-P that are connected by loop range ----------
p <- ggplot(epDF, aes(
        x = cell, y = inLoopPercent, fill = cell, color = cell == "GM12878",
        label = round(inLoopPercent, 2)
      )) + 
  geom_bar(stat = "identity") + 
  geom_text(angle = 90, hjust = 1.1) +
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("black", "red"))
ggsave(paste0(outPrefix, ".E-P_with_loop.percent.batrplot.pdf"), w = 6, h = 3)

# -------------- plot percent of E-P that are connected by extended loops-------
p <- ggplot(epDF, aes(
  x = cell, y = inExtLoopPercent, fill = cell, color = cell == "GM12878",
  label = round(inExtLoopPercent, 2)
)) + 
  geom_bar(stat = "identity") + 
  geom_text(angle = 90, hjust = 1.1) +
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("black", "red"))
ggsave(paste0(outPrefix, ".E-P_with_extLoop.percent.batrplot.pdf"), w = 6, h = 3)

# -------------- plot percent of loops that connect E-P associations -----------
p <- ggplot(epDF, 
            aes(x = cell, y = loopWithEPpercent, fill = cell, 
                color = cell == "GM12878",
                label = round(loopWithEPpercent, 2)
            )) + 
  geom_bar(stat = "identity") + 
  geom_text(angle = 90, hjust = 1.1) +
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("black", "red"))
ggsave(paste0(outPrefix, ".loops_connecting_E-Ps.percent.batrplot.pdf"), w = 6, h = 3)

# ------ plot percent of extended loops that connect E-P associations -----------
p <- ggplot(epDF, 
            aes(x = cell, y = extLoopWithEPpercent, fill = cell, 
                color = cell == "GM12878",
                label = round(extLoopWithEPpercent, 2)
            )) + 
  geom_bar(stat = "identity") + 
  geom_text(angle = 90, hjust = 1.1) +
  theme_bw() + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("black", "red"))
ggsave(paste0(outPrefix, ".extLoops_connecting_E-Ps.percent.batrplot.pdf"), w = 6, h = 3)
