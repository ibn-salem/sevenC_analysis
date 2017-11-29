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
require(ggsignif)     # to calculate p-values on plots
source("R/chromloop.functions.R")
source("R/gr_associations.R")

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE

MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1

outPrefix <- "results/v05_NeuroD1_TSS_DIST"

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
peaksGR_macs <- rtracklayer::import.bed(peak_files[[1]],
               seqinfo = seqInfo_mm9,
               extraCols = extraCols_narrowPeak)

peaksGR_TSS <-  rtracklayer::import.bed(peak_files[[2]], seqinfo = seqInfo_mm9)
peaksGR_DIST <-  rtracklayer::import.bed(peak_files[[3]], seqinfo = seqInfo_mm9)

peakn <- tibble(
  source = names(peak_files),
  n = map_int(list(peaksGR_macs, peaksGR_TSS, peaksGR_DIST), length),
  files = unlist(peak_files),
) %>% 
  write_tsv(paste0(outPrefix, ".n_peaks.tsv"))

peaksGR <- c(peaksGR_TSS, peaksGR_DIST)

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
tssGR$DE_1.5 <- factor(
  tssGR$padj <= 0.05 & (tssGR$log2FoldChange >= 1.5 | tssGR$log2FoldChange <= -1.5),
  c(FALSE, TRUE), c("Not DE", "DE")
  )
tssGR$DE_1 <- factor(
  tssGR$padj <= 0.05 & (tssGR$log2FoldChange >= 1 | tssGR$log2FoldChange <= -1),
  c(FALSE, TRUE), c("Not DE", "DE")
)
tssGR$DE_0.5 <- factor(
  tssGR$padj <= 0.05 & (tssGR$log2FoldChange >= 0.5 | tssGR$log2FoldChange <= -0.5),
  c(FALSE, TRUE), c("Not DE", "DE")
)
tssGR$DE <- tssGR$DE_1.5

write_rds(tssGR, paste0(outPrefix, ".tssGR"))
# tssGR <- read_rds(paste0(outPrefix, ".tssGR"))

#-------------------------------------------------------------------------------
# Get bound gens
#-------------------------------------------------------------------------------


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
isLinked <- function(linkDF, len){
  # 1:len %in% linkDF$gr1
  is.element(seq(len), linkDF$gr1)
}


maxgap_sizes <- c(0, 10^3, 5*10^3, 10^4, 5*10^4, 10^5)

regPredDF <- tibble(
  maxgap = maxgap_sizes,
  ) %>% 
  mutate(
    overlap = map(maxgap, ~ overlapsAny(tssGR, peaksGR, maxgap = .)),
    loop_maxgap = map(
      map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, maxgap = .)),
      isLinked, length(tssGR)
    ),
    loop_inner_maxgap = map(
      map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, inner_maxgap = .)),
      isLinked, length(tssGR)
    ),
    loop_inner_outer = map(
      map(maxgap, ~ linkRegions(tssGR, peaksGR, subGI, inner_maxgap = ., outer_maxgap = .)),
      isLinked, length(tssGR)
    ),
    loop_inLoop = map(
      map(maxgap, ~ linkRegionsInLoops(tssGR, peaksGR, subGI, maxgap = .)), 
      isLinked, length(tssGR)
    )
  ) %>% 
  dplyr::select(maxgap, overlap, starts_with("loop_")) %>% 
  gather(key = method, value = pred, overlap, starts_with("loop_"))

write_rds(regPredDF, paste0(outPrefix, ".regPredDF.rds"))
# regPredDF <- read_rds( paste0(outPrefix, ".regPredDF.rds"))

# calculate contingency tables and performance evaluations
de_threshold <- c("DE_0.5", "DE_1", "DE_1.5")
regGeneList <- map(de_threshold, function(de_col){
  regPredDF %>% 
    mutate(
      contab = map(pred, ~ getContab(mcols(tssGR)[, de_col] == "DE", .x)),
      bc = map(contab, networkBMA::scores)
    ) %>% 
    unnest(bc) %>% 
    dplyr::select(-pred, -contab)
})
names(regGeneList) <- de_threshold


regGeneDF <- bind_rows(regGeneList, .id = "de_threshold")

write_tsv(regGeneDF, paste0(outPrefix, ".regGeneDF.tsv"))
# regGeneDF <- read_tsv(paste0(outPrefix, ".regGeneDF.tsv"))
#-------------------------------------------------------------------------------
# Plot performance
#-------------------------------------------------------------------------------

tidyRegDF <- regGeneDF %>% 
  dplyr::select(-expected, -`O/E`) %>% 
  mutate(maxgap = as.factor(maxgap)) %>% 
  gather(key = "metric", value = "performance", TPR:ACC) %>% 
  filter(metric %in% c("F1", "precision", "sensitivity", "specificity", "MCC"))

p <- ggplot(tidyRegDF, aes(x = maxgap, y = performance, fill = method, label = round(performance, 2))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(hjust = -.1, position = position_dodge(.9), angle = 90) +
  # facet_wrap(~ metric) + #, scales = "free_y"
  facet_grid(metric ~ de_threshold) +
  ylim(0, 1.2) +
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".regGeneDF.binary_classification.barplot.pdf"), w = 12, h = 12)

p <- ggplot(tidyRegDF, aes(x = method, y = performance, fill = method, label = round(performance, 2))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(hjust = -.1, position = position_dodge(.9), angle = 90) +
  # facet_wrap(~ metric) + #, scales = "free_y"
  facet_grid(metric + de_threshold ~ maxgap) +
  ylim(0, 1.2) +
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(5, "Set1"))

ggsave(paste0(outPrefix, ".regGeneDF.binary_classification_by_maxgap.barplot.pdf"), w = 12, h = 12)

# plot ROC like kurve
rocDF <- tidyRegDF %>% 
  filter(metric %in% c("sensitivity", "specificity")) %>% 
  spread(key = metric, value = performance)

p <- ggplot(rocDF, aes(x = 1 - specificity, y = sensitivity, color = method)) +
  geom_abline(intercept = 0, slope = 1, color = "black") + 
  geom_point() +
  geom_line() + 
  facet_grid(de_threshold ~ .) + 
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = brewer.pal(5, "Set1"), 
                     guide = guide_legend(ncol = 2, title.position = "top"))
ggsave(paste0(outPrefix, ".regGeneDF.ROC_by_de_threshold.pdf"), w = 3, h = 9)
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
  facet_grid(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size = 10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = brewer.pal(5, "Set1"))
ggsave(paste0(outPrefix, ".DEvsRegDF.padj_vs_binding.boxplot.pdf"), w = 6, h = 6)


# abs log fold change
p <- ggplot(DEvsRegDF, aes(x = method, y = abs(log2FoldChange), color = pred)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size = 10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0, .5) +
  scale_color_manual(values = brewer.pal(5, "Set1"))
ggsave(paste0(outPrefix, ".DEvsRegDF.absLogFC_vs_binding.boxplot.pdf"), w = 6, h = 6)

# abs log fold change
p <- ggplot(DEvsRegDF, aes(x = pred, y = abs(log2FoldChange), color = pred)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_signif(comparisons = list(c(1, 2)), map_signif_level = TRUE) +
  geom_signif(comparisons = list(c(1, 2)), color = "black", y = 0.8) +
  facet_grid(maxgap ~ method) +
  theme_bw() +
  theme(
    text = element_text(size = 10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0, 1) +
  scale_color_manual(values = brewer.pal(5, "Set1"))
ggsave(paste0(outPrefix, ".DEvsRegDF.absLogFC_by_method.boxplot.pdf"), w = 6, h = 12)


# raw log fold change
p <- ggplot(DEvsRegDF, aes(x = method, y = log2FoldChange, color = pred)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ maxgap) +
  theme_bw() +
  theme(
    text = element_text(size = 10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(-2.5, 3) +
  scale_color_manual(values = brewer.pal(5, "Set1"))
ggsave(paste0(outPrefix, ".DEvsRegDF.LogFC_vs_binding.boxplot.pdf"), w = 6, h = 6)

