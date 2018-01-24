################################################################################
# Analysis of predictd chromatin looping interactions using the sevenC tool
################################################################################


require(sevenC)    # devtools::install_github("ibn-salem/sevenC")
require(rtracklayer)  # to import() BED files
require(EnsDb.Hsapiens.v75)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for human genes
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(feather)      # for efficient storing of data.frames
require(networkBMA)
source("R/sevenC.functions.R")
source("R/gr_associations.R")

# 0) Set parameter --------------------------------------------------------

MIN_MOTIF_SIG <- 6
WINDOW_SIZE <- 1000
BIN_SIZE <- 1

PreviousOutPrefix <- "results/v03_screen_TF_qfraq.motifSig6_w1000_b1"
outPrefix <- "results/v05_GR"

# metadata file
metaFile <- "data/GR_CHIP-SEQ_FILES.tsv"
# peak_files <- c("data/Vockley2016/ChIP-seq/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed",
#                 "data/Vockley2016/ChIP-seq/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed")
peak_file = "data/GR_BAM/ENCFF080RMP.bed"

GR_fc_file = "data/GR_BAM/ENCFF678TFX.bigWig"
# GR_fc_file = "data/ENCODE/Experiments/ENCFF563WQV.bam-qfrags_allChr_chip.bed.sorted.bedGraph.bw"

DE_file = "data/Vockley2016/RNA-seq/countDF.DESeq2.results.tsv"

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

# parse ucscMeta file
meta <- read_tsv(metaFile) %>% 
  filter(output_type == "qfraq")


load(paste0(PreviousOutPrefix, ".gi.tmp.Rdata"))

# reforamt gi data
mcols(gi) <- as.data.frame(mcols(gi)) %>%
  mutate(
    id = 1:nrow(.),
    loop = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop")) 
  ) %>% 
  dplyr::select(id, loop, everything()) %>% 
  DataFrame()


# parse peaks
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
peaksGR <- rtracklayer::import(peak_file, format = "BED",
                        extraCols = extraCols_narrowPeak, seqinfo = seqinfo(gi))

# differntial expression
DE_results <- read_tsv(DE_file, col_types = cols(
  baseMean = col_double(),
  log2FoldChange = col_double(),
  lfcSE = col_double(),
  stat = col_double(),
  pvalue = col_double(),
  padj = col_double()
)) %>% 
  filter(!is.na(ensg))


# load default model parameters
modelDF <- read_tsv(paste0(PreviousOutPrefix, ".bestNModelDF.tsv"), col_types = cols(
  term = col_character(),
  estimate_mean = col_double(),
  estimate_median = col_double(),
  estimate_sd = col_double()
))

# load cutoff value
cutoffDF <- read_tsv(paste0(PreviousOutPrefix, ".topNf1ModelDF.tsv"), col_types = cols(
  mean_max_cutoff = col_double(),
  median_max_cutoff = col_double(),
  sd_max_cutoff = col_double()
))

#-------------------------------------------------------------------------------
# Annotae with coverage and correlation -------------------------------------
#-------------------------------------------------------------------------------


# check if all files exists
meta <- meta %>%
  filter(file.exists(path))


# iterate over all ChIP-seq sata sets
for (i in seq_len(nrow(meta))) {

  message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
  
  gi <- addCor(gi, bwFile = meta$path[i], name = meta$name[i])

}

# # combine replicates
# treatment_names <- c("cor_ChIP_DEX_rep1", "cor_ChIP_DEX_rep2") 
# mcols(gi)[, "cor_ChIP_DEX"] <- mcols(gi)[, treatment_names[1]] + mcols(gi)[, treatment_names[2]] / 2

# # add coverage for GR fold/change
gi <- gi %>% 
  addCor(GR_fc_file, name = "GR_foldChange")

#-------------------------------------------------------------------------------
# Predict interactions
#-------------------------------------------------------------------------------

# write_rds(gi, paste0(outPrefix, "gi_tmp.rds"))


# write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

# for (nameStr in meta$name) {
# for (nameStr in c("GR_foldChange")) {
for (nameStr in c(meta$name, "GR_foldChange")) {
    
  design <- as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", nameStr))
  
  # mcols(gi)[, paste0("pred_", nameStr)] <- pred_logit(mcols(gi), design, modelDF$estimate_mean)
  # mcols(gi)[, paste0("predBinary_", nameStr)] <-  mcols(gi)[, paste0("pred_", nameStr)] >= cutoffDF$mean_max_cutoff

  gi <- gi %>% 
    predLoops(
      formula = design, 
      betas = modelDF$estimate_mean, 
      colname = paste0("pred_", nameStr),
      cutoff = NULL)
  mcols(gi)[, paste0("predBinary_", nameStr)] <-  mcols(gi)[, paste0("pred_", nameStr)] >= cutoffDF$mean_max_cutoff
}



#-------------------------------------------------------------------------------
# save gi with predictions
#-------------------------------------------------------------------------------
write_rds(gi, paste0(outPrefix, "gi_GR_pred.rds"))
# gi <- read_rds(paste0(outPrefix, "gi_GR_pred.rds"))

#===============================================================================
# Associate Loops to differential expressed genes
#===============================================================================

#-------------------------------------------------------------------------------
# get genes
#-------------------------------------------------------------------------------
edb <- EnsDb.Hsapiens.v75
## Change the seqlevels style form Ensembl (default) to UCSC:
seqlevelsStyle(edb) <- "UCSC"
# options(ensembldb.seqnameNotFound=NA)
# options(ensembldb.seqnameNotFound = "ORIGINAL")


# genesGR <- genes(edb, filter = SeqNameFilter(seqnames(seqinfo(gi))), seqinfo = seqinfo(gi))
# genesGR <- genes(edb, filter = ~ seq_name %in% seqnames(seqinfo(gi)))
genesGR <- genes(edb, filter = ~ gene_id %in% DE_results$ensg)

# change seqinfo ojbect
seqlevels(genesGR) <- seqlevels(gi)
seqinfo(genesGR) <- seqinfo(gi)

# filter and reorder DE results
deDF <- DE_results %>% 
  slice(match(genesGR$gene_id, ensg))

# add as annotation to genesGR
mcols(genesGR) <- cbind(mcols(genesGR), deDF)

# define DE genes
genesGR$DE <- factor(
  genesGR$padj <= 0.05 & 
    (genesGR$log2FoldChange >= .5 | genesGR$log2FoldChange <= -.5),
  c(FALSE, TRUE),
  c("Not DE", "DE"))

tssGR <- resize(genesGR, width = 1, fix = "start")

#-------------------------------------------------------------------------------
# Get bound gens
#-------------------------------------------------------------------------------

subGI <- gi[!is.na(gi$predBinary_GR_foldChange) & gi$predBinary_GR_foldChange]

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


regPredDF <- tibble(
  maxgap = c(0, 10^3, 10^4, 10^5, 10^6),
  ) %>% 
  mutate(
    gr = map(maxgap, ~ add_overalp(tssGR, peaksGR, colname = paste0("peak_ovlerap_", .), maxgap = .)),
    ovlerap = map2(gr, maxgap, ~ mcols(.x)[, paste0("peak_ovlerap_", .y)]),
    loop_maxgap = map(map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, maxgap = .)), isLinked, length(tssGR)),
    loop_inner_maxgap = map(map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, inner_maxgap = .)), isLinked, length(tssGR)),
    loop_inner_outer = map(map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, inner_maxgap = ., outer_maxgap = .)), isLinked, length(tssGR)),
    loop_inLoop = map(map(maxgap, ~ linkRegionsInLoops(tssGR, peaksGR, subGI, maxgap = .)), isLinked, length(tssGR)),
  ) %>% 
  dplyr::select(maxgap, ovlerap, starts_with("loop_")) %>% 
  gather(key = method, value = pred, ovlerap, starts_with("loop_"))

write_rds(regPredDF, paste0(outPrefix, ".regPredDF.rds"))

# calculate contingency tables and performance evaluations
regGeneDF <- regPredDF %>% 
  mutate(
    contab = map(pred, ~ getContab(tssGR$DE == "DE", .x)),
    bc = map(contab, networkBMA::scores)
  ) %>% 
  unnest(bc) %>% 
  dplyr::select(-pred, -contab)
  
write_tsv(regGeneDF, paste0(outPrefix, ".regGeneDF.tsv"))

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

ggsave(paste0(outPrefix, ".regGeneDF.binary_classification.barplot.pdf"), w = 12, h = 6)

#-------------------------------------------------------------------------------
# compare predicted regulated genes to log10(p-value) of expression
#-------------------------------------------------------------------------------
DEvsRegDF <- regPredDF %>% 
  mutate(
    padj = list(genesGR$padj),
    log2FoldChange = list(genesGR$log2FoldChange)
    ) %>% 
  unnest(pred, padj, log2FoldChange)

p <- ggplot(DEvsRegDF, aes(x = pred, y = padj, color = method)) +
  geom_boxplot() +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".DEvsRegDF.padj_vs_binding.boxplot.pdf"), w = 6, h = 6)

# abs log fold change
p <- ggplot(DEvsRegDF, aes(x = pred, y = abs(log2FoldChange), color = method)) +
  geom_boxplot() +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".DEvsRegDF.absLogFC_vs_binding.boxplot.pdf"), w = 6, h = 6)

