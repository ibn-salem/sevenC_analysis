################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)   # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)
require(stringr)


# 0) Set parameter --------------------------------------------------------


# work only on subset here:

useTFs <- c(
  "ZNF143",
  "STAT1",
  "ZNF274",
  "ZNF384",
  "CTCF",
  "RAD21",
  "POLR2A",
  "NFYB"
)


# 1) select motifs and parse input data -----------------------------------
ancGR <- chromloop::motif.hg19.CTCF

# parse metatdata table for input files
all_metadata <- read_tsv("data/ENCODE/metadata.flt.tsv")

# filter input data
meta <- all_metadata %>% 
  filter(TF %in% str_to_upper(useTFs)) %>%
  mutate(TF = parse_factor(TF, levels = useTFs)) %>%
  arrange(TF)

# DEBUG: use only the first two ChIP-seq data sets
meta <- head(meta, 2)

# 2) annotae with coverage ------------------------------------------------
covRle <- lapply(meta$filePath, parseBigWigToRle, seqInfo = seqinfo(ancGR))

for (i in seq_along(meta$TF)) {
  ancGR <- chromloop::addCovToGR(ancGR, covRle[[i]], colname = paste0("cov_", meta$TF[i]))
}

# 3) build InteractonSet object -------------------------------------------
gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)

# 4) compute correlations -------------------------------------------------
gi <- chromloop::applyToCloseGI(gi, cor, str_c("cor_", meta$TF[1]))

# 5) add ture loops -------------------------------------------------------

# 6) predict loops --------------------------------------------------------

# 7)  Analyse performace --------------------------------------------------

