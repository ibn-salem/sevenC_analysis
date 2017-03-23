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
meta <- read_tsv("data/ENCODE/metadata.flt.tsv")

# # filter input data
# meta <- meta %>% 
#   filter(TF %in% useTFs) %>%
#   mutate(TF = parse_factor(TF, levels = useTFs)) %>%
#   arrange(TF)

# filter input data
meta <- meta %>% 
  filter(TF %in% useTFs) %>%
  filter(`Output type` == "fold change over control") %>%
  mutate(TF = parse_factor(TF, intersect(useTFs, .$TF))) %>%
  arrange(TF)

# # DEBUG use only two examples here:
# meta <- meta %>%
#   filter(TF %in% c("POLR2A", "CTCF"))

# 2) annotae with coverage ------------------------------------------------

# covRle <- lapply(meta$filePath, parseBigWigToRle, seqInfo = seqinfo(ancGR))

for (i in seq_along(meta$TF)) {
  ancGR <- chromloop::addCovToGR(ancGR, meta$filePath[i], window = 10, 
                                 colname = paste0("cov_", meta$TF[i]))
}

# 3) build InteractonSet object -------------------------------------------
gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)

# 4) compute correlations -------------------------------------------------
for (i in seq_along(meta$TF)) {
  gi <- chromloop::applyToCloseGI(gi, 
                                  datcol = paste0("cov_", meta$TF[i]),
                                  fun = cor, 
                                  colname = paste0("cor_", meta$TF[i]))
}

# 5) add ture loops -------------------------------------------------------

# 6) predict loops --------------------------------------------------------

# 7)  Analyse performace --------------------------------------------------

