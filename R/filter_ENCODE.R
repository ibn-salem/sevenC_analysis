

# a script to filter ENCODE data sets by metadta table. 

require(tidyverse)
require(stringr)

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

#-------------------------------------------------------------------------------
# parse meta data and filter
#-------------------------------------------------------------------------------

metaDataFile <- "data/ENCODE/metadata.tsv"

# parse metadata
meta <- read_tsv(metaDataFile)

# add TF name
df <- meta %>%
  mutate(TF = str_replace(`Experiment target`, "-human", "")) %>%
  select(`File accession`, TF, `Output type`, 
         `Biological replicate(s)`, Lab, everything())

# utils:::format.object_size(x, format="auto")

plotDF <- df %>%
  filter(`Biosample term name` == "GM12878") %>%
  filter(Lab == "ENCODE Processing Pipeline") %>%
  mutate(size = parse_double(Size)) %>%
  group_by(`Output type`, `Biological replicate(s)`) %>%
  summarise(
      mean_size = mean(size, na.rm=TRUE),
      n = n()
  )

# utils:::format.object_size(plotDF$mean_size, "Mb")

plotDF %>%
  ggplot(aes(x=`Output type`, y=n, fill=`Biological replicate(s)`)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))


plotDF %>%
  ggplot(aes(x=`Output type`, y=mean_size, fill=`Biological replicate(s)`)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

  
# filter for "signal" as output with combined replicate (is NA here) and only GM12878 cell line

flt <- df %>% 
  filter(`Output type` %in% c("raw signal", "fold change over control", "signal")) %>%
  filter(`Output type` %in% c("raw signal", "fold change over control")) %>%
  filter(`Biological replicate(s)` == "1" | is.na(`Biological replicate(s)`)) %>% 
  filter(`Biosample term name` == "GM12878") %>%
  filter(TF %in% useTFs)
  

# is.na(`Biological replicate(s)`) | 
# %>%
#   distinct(TF, Lab, .keep_all = TRUE)



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


