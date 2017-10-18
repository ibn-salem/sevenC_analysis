################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
# by predicting regulated genes as validated by TF perturbation experiments
################################################################################


library(chromloop)    # devtools::install_github("ibn-salem/chromloop")
library(rtracklayer)  # to import() BED files
library(EnsDb.Hsapiens.v75)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for human genes
library(networkBMA)
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(modelr)       # for tidy modeling
library(precrec)      # for ROC and PRC curves
library(RColorBrewer)   # for nice colors
library(rsample)
library(pryr) # for object_size()
library(feather)      # for efficient storing of data.frames
library(multidplyr)   # for partition() and collect() to work in parallel
library(readxl)       # to read excel files
source("R/chromloop.functions.R")
source("R/gr_associations.R")


# 0) Set parameter --------------------------------------------------------

# MIN_MOTIF_SIG <- 6
# WINDOW_SIZE <- 1000
# BIN_SIZE <- 1

# PreviousOutPrefix <- "results/v03_screen_TF_qfraq.motifSig6_w1000_b1"
outPrefix <- "results/v04.TF_perturbation"

# metadata file
# metaFile <- "data/GR_CHIP-SEQ_FILES.tsv"

sample_names <- c("MYC_HeLa", "MYC_MCF7")

peak_files <- paste0("data/Wu2013/", c("ENCFF002CSE.bed", "ENCFF002DBF.bed"))

# peak_file = "data/GR_BAM/ENCFF080RMP.bed"

bigWig_file = paste0("data/Wu2013/", c("ENCFF000XCK.bigWig", "ENCFF000RYM.bigWig"))
# GR_fc_file = "data/ENCODE/Experiments/ENCFF563WQV.bam-qfrags_allChr_chip.bed.sorted.bedGraph.bw"

GS_genes_file = "data/Wu2013/12859_2012_5924_MOESM1_ESM.xls"


# ----------------------- Parse peaks ------------------------------------------
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peaks <- map(peak_files, rtracklayer::import, format = "BED",
                        extraCols = extraCols_narrowPeak, seqinfo = seqinfo(gi))

# ----------------------- Parse gold standard genes from perturbation data -----

sheets_lables = c("GS_MYC_HelaS3", "GS_MYC_MCF7")

GS_genes <- map(sheets_lables, ~ read_excel(GS_genes_file, sheet = .x)) %>% 
               map(pull, Original)
# ----------- Predict Loops for input ChiP-seq data ----------------------------
gi <- prepareCandidates(motif.hg19.CTCF, scoreColname = "sig")

for (i in seq_along(sample_names)) {
  
  nameStr <- sample_names[[i]]
  
  gi <- gi %>% 
    addCor(bwFile = bigWig_file[[i]], name = nameStr)
  
  # define model based on specific TF ChIP-seq data set
  model <- as.formula(paste0("~ dist + strandOrientation + score_min + ", nameStr))
  
  # add preiction score
  gi <- gi %>% 
    predLoops(formula = model, cutoff = NULL, colname	= paste0("pred_", nameStr))
  
  # add binary prediction
  mcols(gi)[, paste0("predBinary_", nameStr)] <-  mcols(gi)[, paste0("pred_", nameStr)] >= chromloop::cutoffBest10
}

write_rds(gi, paste0(outPrefix, "gi_pred.rds"))

#=== Associate Loops to differential expressed genes ===========================

#--- get genes ----------------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genesGR <- genes(txdb)

# check if all gold standard genes are contained in txdb
map(GS_genes, ~ mean(.x %in% genesGR$gene_id))

# annotate genesGR with DE info
for (i in seq_along(GS_genes)) {
  
  de <- genesGR$gene_id %in% GS_genes[[i]]
  deFact <- factor(de, c(FALSE, TRUE), c("Not DE", "DE"))
  
  mcols(genesGR)[ ,paste0("DE_", sample_names[[i]])] <- deFact
}

tssGR <- resize(genesGR, width = 1, fix = "start")

# --------------- Some functions -----------------------------------------------
#' Add colum to GenomicRanges object indicating any overlap
#' 
add_overalp <- function(tssGR, peaksGR, colname = "peak_ovwerlap", ...){
  
  mcols(tssGR)[, colname] <- GenomicRanges::countOverlaps(
    tssGR, peaksGR, ignore.strand = FALSE, ...) > 0
  return(tssGR)
  
}


#' Get contingency table as tidy tibble from.
#'
#' @param true Logical vector of true labels (goldstandard)
#' @param pred Logical vector of predicted lables
#' @return A tibble with the number of true postivies (TP), false postivives
#'   (FP), false negatives (FN), and true negatives (TN) in columns.
getContab <- function(true, pred){
  
  tab <- count(
    tibble(
      true = true,
      pred = pred
    ), 
    true, pred
  ) %>% 
    arrange(true, pred)
  
  tibble(
    TN = tab$n[[1]],
    FP = tab$n[[2]],
    FN = tab$n[[3]],
    TP = tab$n[[4]]
  )
}

# helper function
isLinked <- function(linkDF, len){1:len %in% linkDF$gr1}


# --------------- Get bound gens -----------------------------------------------

for (i in seq_along(sample_names)) {
  
  nameStr <- sample_names[[i]]
  
  binaryPred <- mcols(gi)[, paste0("predBinary_", nameStr)]
  
  subGI <- gi[!is.na(binaryPred) & binaryPred]

  regPredDF <- tibble(
    maxgap = c(0, 10^3, 10^4, 10^5, 10^6),
  ) %>% 
    mutate(
      gr = map(maxgap, ~ add_overalp(tssGR, peaks[[i]], colname = paste0("peak_ovlerap_", .), maxgap = .)),
      ovlerap = map2(gr, maxgap, ~ mcols(.x)[, paste0("peak_ovlerap_", .y)]),
      loop_maxgap = map(map(maxgap, ~ linkRegions(tssGR, peaks[[i]], subGI, maxgap = .)), isLinked, length(tssGR)),
      loop_inner_maxgap = map(map(maxgap, ~ linkRegions(tssGR, peaks[[i]], subGI, inner_maxgap = .)), isLinked, length(tssGR)),
      loop_inner_outer = map(map(maxgap, ~ linkRegions(tssGR, peaks[[i]], subGI, inner_maxgap = ., outer_maxgap = .)), isLinked, length(tssGR)),
      loop_inLoop = map(map(maxgap, ~ linkRegionsInLoops(tssGR, peaks[[i]], subGI, maxgap = .)), isLinked, length(tssGR)),
    ) %>% 
    dplyr::select(maxgap, ovlerap, starts_with("loop_")) %>% 
    gather(key = method, value = pred, ovlerap, starts_with("loop_"))
  
  write_rds(regPredDF, paste0(outPrefix, ".", nameStr, ".regPredDF.rds"))
  
  # calculate contingency tables and performance evaluations
  deGenes <- mcols(tssGR)[, paste0("DE_", nameStr)] == "DE"
  
  regGeneDF <- regPredDF %>% 
    mutate(
      contab = map(pred, ~ getContab(deGenes, .x)),
      bc = map(contab, networkBMA::scores)
    ) %>% 
    unnest(bc) %>% 
    dplyr::select(-pred, -contab)
  
  write_tsv(regGeneDF, paste0(outPrefix, ".", nameStr, ".regGeneDF.tsv"))

  
  #-------------------------------------------------------------------------------
  # Plot performance
  #-------------------------------------------------------------------------------
  
  tidyRegDF <- regGeneDF %>% 
    dplyr::select(-expected, -`O/E`) %>% 
    mutate(maxgap = as.factor(maxgap)) %>% 
    gather(key = "metric", value = "performance", TPR:ACC) %>% 
    filter(metric %in% c("F1", "precision", "recall", "sensitivity", "specificity", "MCC"))
  
  p <- ggplot(tidyRegDF, aes(x = maxgap, y = performance, fill = method, label = round(performance, 2))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(hjust = -.1, position = position_dodge(.9), angle = 90) +
    facet_wrap(~ metric) + #, scales = "free_y"
    ylim(0, 1.2) +
    theme_bw() +
    theme(
      text = element_text(size=10), 
      legend.position = "bottom",
      axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(values = brewer.pal(5, "Set1"))
  
  ggsave(paste0(outPrefix, ".", nameStr, ".regGeneDF.binary_classification.barplot.pdf"), w = 12, h = 6)
  
  
}
  

# #-------------------------------------------------------------------------------
# # compare predicted regulated genes to log10(p-value) of expression
# #-------------------------------------------------------------------------------
# DEvsRegDF <- regPredDF %>% 
#   mutate(padj = list(genesGR$padj)) %>% 
#   unnest(pred, padj)
# 
# p <- ggplot(DEvsRegDF, aes(x = pred, y = padj, color = method)) +
#   geom_boxplot() +
#   facet_wrap(~ maxgap) +
#   theme_bw() +
#   theme(
#     text = element_text(size=10), 
#     legend.position = "bottom",
#     axis.text.x = element_text(angle = 60, hjust = 1)) +
#   scale_color_manual(values = brewer.pal(5, "Set1"))
# 
# ggsave(paste0(outPrefix, ".DEvsRegDF.padj_vs_binding.boxplot.pdf"), w = 6, h = 6)


