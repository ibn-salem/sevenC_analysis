

# a script to filter ENCODE data sets by metadta table. 

require(tidyverse)
require(stringr)

outPrefix <- file.path("results", "data.ENCODE")
dir.create(dirname(outPrefix), showWarnings = FALSE)

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
reportFile <- "data/ENCODE/report.tsv"


# parse report.tsv with metadata for each experiment ----------------------

# parse report containing metadata for each experiment
report <- read_tsv(
  reportFile,
  col_types = cols(
    `Biological replicate` = col_character(),
    `Technical replicate` = col_character()
  ))

report <- report %>%
  mutate(exp_nrep = map_int(str_split(`Biological replicate`, ","), length))

# parse metadata.tsv for each file ----------------------------------------

# parse metadata for each file
meta <- read_tsv(
  metaDataFile,
  col_types = cols(
    `Biological replicate(s)` = col_character(),
    `Technical replicate` = col_character()
  ))

# filter for hg19
meta <- meta %>%
  filter(Assembly == "hg19")

# add TF name from report table
meta <- meta %>%
  left_join(
    report,
    select(report, Accession, `Target label`),
    by = c("Experiment accession" = "Accession")
  ) %>%
  mutate(TF = `Target label`, Lab = Lab.x) %>%
  mutate(rep = map(str_split(`Biological replicate(s)`, ","), parse_integer)) %>%
  mutate(file_nrep = map_int(rep, length)) %>%
  mutate(file_nrep = ifelse(is.na(rep), "Unknown", file_nrep)) %>%
  mutate(size = Size / 1024^2)


# reorder columns to have important column at beginning
df <- meta %>%
  select(`File accession`, TF, `Output type`, 
         rep, file_nrep, exp_nrep, Lab, size, everything())


# utils:::format.object_size(x, format="auto")
plotDF <- df %>%
  filter(`Biosample term name` == "GM12878") %>%
  # filter(Lab == "ENCODE Processing Pipeline") %>%
  group_by(`Output type`, file_nrep) %>%
  summarise(
      mean_size = mean(size, na.rm = TRUE),
      n = n()
  )

#plotDF %>% count(`Output type`)
# utils:::format.object_size(plotDF$mean_size, "Mb")

p <- plotDF %>%
  ggplot(aes(x = `Output type`, y = n, fill = factor(file_nrep))) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_c(outPrefix, ".n_by_outputType_and_nrep.barplot.pdf"))


p <- plotDF %>%
  ggplot(aes(x=`Output type`, y=mean_size, fill = factor(file_nrep))) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(str_c(outPrefix, ".mean_size_by_outputType_and_nrep.barplot.pdf"))

p <- df %>%
  filter(`Biosample term name` == "GM12878") %>%
  ggplot(aes(x=`Output type`, y=size, fill=factor(file_nrep))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(str_c(outPrefix, ".size_by_outputType_and_nrep.boxplot.pdf"))



# Filter for used data sets -----------------------------------------------


  
# filter for "signal" as output with combined replicate (is NA here) and only GM12878 cell line

flt <- df %>% 
  filter(`Output type` %in% c("raw signal", "fold change over control")) %>%
  # filter(`Biological replicate(s)` == "1" | is.na(rep)) %>% 
  # filter(map_lgl(rep, setequal, 1:2) | is.na(rep)) %>% 
  filter(file_nrep == 2 | `Output type` == "raw signal") %>% 
  filter(`Biosample term name` == "GM12878") %>%
  distinct(`Output type`, TF, .keep_all = TRUE) %>%
  mutate(usedTF = TF %in% useTFs) %>%
  arrange(desc(usedTF), desc(`Output type`)) %>%
  mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`)))
  
# get size of files
sizeDF <- flt %>%
  group_by(`Output type`) %>%
  summarize(
    n = n(),
    total_size = sum(size) / 1024,
    mean_size = mean(size)
  )

usedFlt <- flt %>%
  filter(usedTF)


# %>%
#   filter(TF %in% useTFs)
  

# is.na(`Biological replicate(s)`) | 
# %>%
#   distinct(TF, Lab, .keep_all = TRUE)

#-------------------------------------------------------------------------------
# Output filtered URL list and metadata table
#-------------------------------------------------------------------------------

flt %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.flt.txt"),
    col_names = FALSE)

flt %>%
  select(-rep) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.flt.tsv"))

# do the same for filtered files

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


