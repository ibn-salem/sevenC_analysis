

# a script to filter ENCODE data sets by metadta table. 

require(tidyverse)
require(stringr)

#-------------------------------------------------------------------------------
# parse meta data and filter
#-------------------------------------------------------------------------------

metaDataFile <- "data/ENCODE/metadata.tsv"

# parse metadata
meta <- read_tsv(metaDataFile)

# add TF name
df <- meta %>%
  mutate(TF = str_replace(`Experiment target`, "-human", ""))


# filter for "signal" as output with combined replicate (is NA here) and only GM12878 cell line

flt <- df %>% 
  filter(`Output type` == "signal") %>%
  filter(is.na(`Biological replicate(s)`)) %>% 
  filter(`Biosample term name` == "GM12878") %>%
  distinct(TF, .keep_all = TRUE)

# add location
flt$filePath <- file.path("data", "ENCODE", "Experiments", basename(flt$`File download URL`))

#-------------------------------------------------------------------------------
# Output filtered URL list and metadata table
#-------------------------------------------------------------------------------

flt %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.flt.txt"),
    col_names = FALSE)

write_tsv(flt, 
          path = file.path("data", "ENCODE", "metadata.flt.tsv"))

#-------------------------------------------------------------------------------
# download files
#-------------------------------------------------------------------------------
# dir.create(file.path("data", "ENCODE", "Experiments"), showWarnings = FALSE)
# 
# # iterate over files and download each
# # for (i in seq_along(flt$filePath)) {
# for (i in 1:2) {
#     
#   message("INFO: Downloading from: ", flt$`File download URL`[i])
#   
#   download.file(flt$`File download URL`[i], flt$filePath[i], mode = "wb")  
# 
# }


