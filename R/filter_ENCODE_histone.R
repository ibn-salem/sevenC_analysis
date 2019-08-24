

# a script to filter ENCODE data sets by metadta table. 

require(tidyverse)
require(stringr)
require(lubridate)

outPrefix <- file.path("results", "data.ENCODE_hisone")
dir.create(dirname(outPrefix), showWarnings = FALSE)

# Top 5 frequent marks (across all tissues)
# H3K4me3 497
# H3K4me1 414
# H3K36me3 408
# H3K27ac 386
# H3K27me3 379
useMarks <- c(
  "H3K4me3",
  "H3K4me1",
  "H3K36me3",
  "H3K27ac",
  "H3K27me3"
)

#-------------------------------------------------------------------------------
# parse meta data and filter
#-------------------------------------------------------------------------------
metaDataFile <- "data/ENCODE/histone/metadata.tsv"
# reportFile <- "data/ENCODE/report.tsv"

# parse metadata.tsv for each file ----------------------------------------

# parse metadata for each file
meta_raw <- read_tsv(
  metaDataFile,
  col_types = cols(
    `Biological replicate(s)` = col_character(),
    `Technical replicate` = col_character()
  ))

# filter for files of interest
meta <- meta_raw %>% 
  filter(
    `Output type` == "fold change over control",
    `Biological replicate(s)` == "1, 2",
    Assembly == "hg19"
    ) %>% 
  mutate(filePath = file.path("data", "ENCODE", "histone", basename(`File download URL`)))

#-------------------------------------------------------------------------------
# Filter for fold change and all TFs -------------------------------------------------
#-------------------------------------------------------------------------------

write_tsv(meta, file.path("data", "ENCODE", "metadata.meta.tsv"))

meta %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "histone", "URLs.txt"),
    col_names = FALSE)

