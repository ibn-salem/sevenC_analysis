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
peak_files <- c("Vockley2016/ChIP-seq/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed",
                "Vockley2016/ChIP-seq/GSM2095218_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed")

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

# parse ucscMeta file
meta <- read_tsv(metaFile) %>% 
  filter(output_type == "shifted_reads")


load(paste0(PreviousOutPrefix, ".gi.tmp.Rdata"))

# parse peaks
peaks <- map(peak_files, rtracklayer::import.bed, seqinfo = seqinfo(gi))

# load default model parameters
modelDF <- read_tsv(paste0(PreviousOutPrefix, "bestNModelDF.tsv"), col_types = cols(
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


#-------------------------------------------------------------------------------
# Predict interactions
#-------------------------------------------------------------------------------

mcols(gi) <- as.data.frame(mcols(gi)) %>%
  mutate(
    id = 1:nrow(.),
    loop = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop")) 
    ) %>% 
    select(id, loop, everything()) %>% 
  DataFrame()


# write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

for (nameStr in meta$name) {
  
  design <- as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", nameStr))
  
  mcols(gi)[, paste0("pred_", name)] <- pred_logit(mcols(gi), design, modelDF$estimate_mean)
  mcols(gi)[, paste0("predBinary_", name)] <-  mcols(gi)[, paste0("pred_", name)] >= cutoffDF$mean_max_cutoff
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#===============================================================================
#===============================================================================
