################################################################################
# Analysis NeuroD1 ChIP-seq data from Pataskar et al. 2016 to predict loops
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(rtracklayer)  # for import() of BED files
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(EnsDb.Mmusculus.v79)  # for mouse genes
require(org.Mm.eg.db) # for mouse gene ID mapping
require(biomaRt)
require(tidyverse)    # for tidy data
require(readxl)       # for reading .xlsx files
require(stringr)      # for string functions
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(feather)      # for efficient storing of data.frames
require(networkBMA)
source("R/chromloop.functions.R")
source("R/gr_associations.R")

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE

MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1

outPrefix <- "results/v05_NeuroD1"

# topN_models_file = "results/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1.bestNModelDF.tsv"
models_prefix = "results/v05_screen_TF_lfc.motifPval2.5e-06_w1000_b1"

CTCF_motifs_mm9 = "data/JASPAR2018/MA0139.1.mm9.sorted.bed"

peak_files = list(
  MACS = "data/Pataskar2016/NeuroD1.sam.bam_sort.bam_peaks.encodePeak",
  TSS = "data/Pataskar2016/TSS_order.bed",
  DIST = "data/Pataskar2016/DIST_order.bed"
)
  
ChIP_input_file = "data/Pataskar2016/GSE65072_Input.wig"
ChIP_NeuoD1_file = "data/Pataskar2016/GSE65072_NeuroD1_mergedReplicates.wig"

DE_file = "data/Pataskar2016/GSE65072_DESeq_plusDOX_minusDOX.xlsx"

#*******************************************************************************
# Parse CTCF motifs for mouse and build pairs  ---------------------------------
#*******************************************************************************

# get seqInfo
seqInfo_mm9 <- seqinfo(TxDb.Mmusculus.UCSC.mm9.knownGene)
seqInfo_mm10 <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)

# parse motifs
motifDF <- read_tsv(CTCF_motifs_mm9, col_names = FALSE) %>% 
  # revers transfomation -log_10(p) * 100
  mutate(score = X5 / 100 ) %>% 
  filter(score > -log10(MOTIF_PVAL))

motifGR <- GRanges(
    motifDF$X1,
    IRanges(motifDF$X2, motifDF$X3),
    strand = motifDF$X6,
    score = motifDF$score,
    seqinfo = seqInfo_mm9
  ) %>% 
  sort()

# filter out-of-chrom regiosn
out_idx <- GenomicRanges:::get_out_of_bound_index(motifGR)
motifGR <- motifGR[-out_idx]

# prepare coandidates
gi <- prepareCisPairs(motifGR)

#*******************************************************************************
# Parse ChIP-seq signals
#*******************************************************************************
pseudo_count <- 1
regions(gi) <- regions(gi) %>% 
  addCovToGR(ChIP_input_file, colname = "input") %>% 
  addCovToGR(ChIP_NeuoD1_file, colname = "NeuroD1_raw")

input <- mcols(regions(gi))[, "input"] %>% as.list()
chip <- mcols(regions(gi))[, "NeuroD1_raw"] %>% as.list()

fe <- map2(chip, input, ~ (.x + pseudo_count) / (.y + pseudo_count))
mcols(regions(gi))[, "NeuroD1"] <- NumericList(fe)

gi <- addCovCor(gi, datacol = "NeuroD1", colname = "cor_NeuroD1")

# save file for faster reload
write_rds(gi, paste0(outPrefix, ".gi.rds"))
# gi <- read_rds( paste0(outPrefix, ".gi.rds") )

#*******************************************************************************
# Predict loops
#*******************************************************************************
# parse model parameter
model <- read_tsv(paste0(models_prefix, ".bestNModelDF.tsv"))
cutoffDF <- read_tsv(paste0(models_prefix, ".topNf1ModelDF.tsv"))

# predict loops
gi <- gi %>% 
  predLoops(
    formula = ~ dist + strandOrientation + score_min + cor_NeuroD1,
    betas = model$estimate_mean,
    colname = "pred",
    cutoff = NULL
  )

subGI <- gi[!is.na(gi$pred) & gi$pred >= cutoffDF$mean_max_cutoff]


#*******************************************************************************
# parse peaks
#*******************************************************************************
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
peaksGR <- rtracklayer::import.bed(peak_files[[1]],
               seqinfo = seqInfo_mm9,
               extraCols = extraCols_narrowPeak)


#*******************************************************************************
# differntial expression
#*******************************************************************************
DE_results <- read_excel(DE_file) %>% 
  rename(GENE_NAME = `GENE NAME`) %>% 
  mutate(
    foldChange = parse_double(foldChange),
    log2FoldChange = parse_double(log2FoldChange),
    pval = parse_double(pval),
    padj = parse_double(padj)
  )

#===============================================================================
# Associate Loops to differential expressed genes
#===============================================================================
mart <- useMart("ENSEMBL_MART_ENSEMBL", 
                dataset = "mmusculus_gene_ensembl", 
                host = "http://may2012.archive.ensembl.org")

# filters = "external_gene_name", 
# values = symbols,

genesDF <- getBM(
  mart = mart,
  attributes = c("external_gene_id", "ensembl_gene_id", "chromosome_name", 
                 "strand", "transcript_start", "transcript_end")
  ) %>% 
  as.tibble()

tssDF <- genesDF %>%
  mutate(transcript_length = transcript_end - transcript_start) %>% 
  arrange(external_gene_id, ensembl_gene_id, desc(transcript_length)) %>% 
  distinct(external_gene_id, ensembl_gene_id, .keep_all = TRUE)

deDF <- DE_results %>% 
  inner_join(tssDF, by = c("GENE_NAME" = "external_gene_id")) %>% 
  filter(paste0("chr", chromosome_name) %in% seqnames(seqInfo_mm9))


tssGR <- GRanges(
  paste0("chr", deDF$chromosome_name),
  IRanges(deDF$transcript_start, deDF$transcript_end),
  strand = ifelse(deDF$strand == 1, "+", "-"),
  seqinfo = seqInfo_mm9
)

tssGR <- resize(tssGR, width = 1, fix = "start")
mcols(tssGR) = DataFrame(deDF)


# define DE genes
tssGR$DE <- factor(
  tssGR$padj <= 0.05 & 
    (tssGR$log2FoldChange >= 1.5 | tssGR$log2FoldChange <= -1.5),
  c(FALSE, TRUE),
  c("Not DE", "DE"))

write_rds(tssGR, paste0(outPrefix, ".tssGR"))

#-------------------------------------------------------------------------------
# Get bound gens
#-------------------------------------------------------------------------------

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
# regPredDF <- read_rds( paste0(outPrefix, ".regPredDF.rds"))

# calculate contingency tables and performance evaluations
regGeneDF <- regPredDF %>% 
  mutate(
    contab = map(pred, ~ getContab(tssGR$DE == "DE", .x)),
    bc = map(contab, networkBMA::scores)
  ) %>% 
  unnest(bc) %>% 
  dplyr::select(-pred, -contab)
  
write_tsv(regGeneDF, paste0(outPrefix, ".regGeneDF.tsv"))
# regGeneDF <- read_tsv(paste0(outPrefix, ".regGeneDF.tsv"))
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
    padj = list(tssGR$padj),
    log2FoldChange = list(tssGR$log2FoldChange)
    ) %>% 
  unnest(pred, padj, log2FoldChange)


p <- ggplot(DEvsRegDF, aes(x = pred, y = padj, color = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size = 10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".DEvsRegDF.padj_vs_binding.boxplot.pdf"), w = 6, h = 6)


# abs log fold change
p <- ggplot(DEvsRegDF, aes(x = pred, y = abs(log2FoldChange), color = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".DEvsRegDF.absLogFC_vs_binding.boxplot.pdf"), w = 6, h = 6)

