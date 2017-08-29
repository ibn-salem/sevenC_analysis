################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for human genes
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(rtracklayer)  # to import() BED files
require(rsample)
require(pryr) # for object_size()
require(feather)      # for efficient storing of data.frames
require(multidplyr)   # for partition() and collect() to work in parallel
source("R/chromloop.functions.R")


# 0) Set parameter --------------------------------------------------------

MIN_MOTIF_SIG <- 6
WINDOW_SIZE <- 1000
BIN_SIZE <- 1

PreviousOutPrefix <- "results/v03_screen_TF_qfraq.motifSig6_w1000_b1"
outPrefix <- "results/v03_screen_TF_qfraq.motifSig6_w1000_b1.GR"

# metadata file
metaFile <- "data/GR_CHIP-SEQ_FILES.tsv"
peak_files <- c("data/Vockley2016/ChIP-seq/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed",
                "data/Vockley2016/ChIP-seq/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed")

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
peaks <- map(peak_files, rtracklayer::import.bed, seqinfo = seqinfo(gi))

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
  
  # add coverage  
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

# add coverage for GR fold/change
regions(gi) <- chromloop::addCovToGR(regions(gi), GR_fc_file, 
  window = WINDOW_SIZE,
  bin_size = BIN_SIZE,
  colname = "cov_GR_foldChange"
)

# add correlations
gi <- gi %>% 
  chromloop::applyToCloseGI(
    datcol = "cov_GR_foldChange", 
    fun = cor, colname = "cor_GR_foldChange")


#-------------------------------------------------------------------------------
# Predict interactions
#-------------------------------------------------------------------------------

# write_rds(gi, paste0(outPrefix, "gi_tmp.rds"))


# write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

# for (nameStr in meta$name) {
for (nameStr in c("GR_foldChange")) {
    
  design <- as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", nameStr))
  
  mcols(gi)[, paste0("pred_", nameStr)] <- pred_logit(mcols(gi), design, modelDF$estimate_mean)
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
require(EnsDb.Hsapiens.v75)
require(networkBMA)

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

peaksGR <- peaks[[1]]

#' Add colum to GenomicRanges object indicating any overlap
#' 
add_overalp <- function(tssGR, peakGR, colname = "peak_ovwerlap", ...){
  
  mcols(tssGR)[, colname] <- GenomicRanges::countOverlaps(
    tssGR, peaksGR, ignore.strand = FALSE, ...) > 0
  return(tssGR)

  }

tssGR <- add_overalp(tssGR, peaksGR, colname = "peak_ovwerlap", maxgap = 0)

linksDF <- chromloop::linkRegions(tssGR, peaksGR, subGI) %>% 
  as_tibble()

tssGR$peak_loop <- seq(length(tssGR)) %in% linksDF$gr1

tab <- table(tssGR$DE, tssGR$peak_loop)

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


regPredDF <- tibble(
  maxgap = c(0, 10^3, 10^4, 10^5, 10^6),
) %>% 
  mutate(
    gr = map(maxgap, ~ add_overalp(tssGR, peakGR, colname = paste0("peak_ovlerap_", .), maxgap = .)),
    peak_ovlerap = map2(gr, maxgap, ~ mcols(.x)[, paste0("peak_ovlerap_", .y)]),
    linkDF = map(maxgap, ~ chromloop::linkRegions(tssGR, peaksGR, subGI, maxgap = .)),
    peak_loop = map(linkDF, ~ seq(length(tssGR)) %in% .$gr1),
  ) %>% 
  select(maxgap, peak_ovlerap, peak_loop) %>% 
  gather(key = method, value = pred, peak_ovlerap, peak_loop)

regGeneDF <- regPredDF %>% 
  mutate(
    contab = map(pred, ~ getContab(tssGR$DE == "DE", .x)),
    bc = map(contab, networkBMA::scores)
  ) %>% 
  unnest(bc) %>% 
  select(-pred, -contab)
  


write_tsv(regGeneDF, paste0(outPrefix, ".regGeneDF.tsv"))

#-------------------------------------------------------------------------------
# Plot performance
#-------------------------------------------------------------------------------

tidyRegDF <- regGeneDF %>% 
  select(-expected, -`O/E`) %>% 
  mutate(maxgap = as.factor(maxgap)) %>% 
  gather(key = "metric", value = "performance", TPR:ACC) %>% 
  filter(metric %in% c("F1", "precision", "recall", "sensitivity", "specificity", "MCC"))

p <- ggplot(tidyRegDF, aes(x = method, y = performance, fill = maxgap, label = round(performance, 2))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(hjust = "inward", position = position_dodge(.9), angle = 90) +
  facet_wrap(~ metric) + 
  theme_bw() +
  theme(
    text = element_text(size=10), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave(paste0(outPrefix, ".regGeneDF.binary_classification.barplot.pdf"), w = 6, h = 6)


