

# a script to filter ENCODE data sets by metadta table. 

require(tidyverse)
require(stringr)
require(lubridate)

outPrefix <- file.path("results", "data.ENCODE")
dir.create(dirname(outPrefix), showWarnings = FALSE)

useTFs <- c(
  "RAD21",
  "CTCF",
  "ZNF143",
  "STAT1",
  "EP300",
  "POLR2A"
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
  skip = 1,
  col_types = cols(
    `Biological replicate` = col_character(),
    `Technical replicate` = col_character()
  ))

report <- report %>%
  mutate(exp_nrep = map_int(str_split(`Biological replicate`, ","), length))

# parse metadata.tsv for each file ----------------------------------------

# parse metadata for each file
meta_raw <- read_tsv(
  metaDataFile,
  col_types = cols(
    `Biological replicate(s)` = col_character(),
    `Technical replicate` = col_character()
  ))

# add TF name from report table, 
# add number of replicates
# normalize size in Mb
# and seleect importat colums to beginning
meta <- meta_raw %>%
  left_join(
    report,
    # select(report, Accession, `Target label`),
    by = c("Experiment accession" = "Accession")
  ) %>%
  mutate(TF = `Target label`, Lab = Lab.x) %>%
  # mutate(rep = map(str_split(`Biological replicate(s)`, ","), parse_integer)) %>%
  mutate(rep = map(str_split(`Biological replicate(s)`, ","), parse_integer)) %>%
  mutate(file_nrep = map_int(rep, length)) %>%
  mutate(file_nrep = ifelse(is.na(rep), "Unknown", file_nrep)) %>%
  mutate(size = Size / 1024^2) %>% 
  mutate(date = ymd(`Experiment date released`)) %>% 
  select(`File accession`, TF, `Output type`, 
         rep, file_nrep, exp_nrep, `File format`, Lab, size, date, everything()) 


# filter for bigWig
# filter for GM12878
# filter for hg19
df <- meta %>%
  filter(`File format` == "bigWig") %>% 
  filter(Assembly == "hg19") %>% 
  filter(`Biosample term name` == "GM12878")



# utils:::format.object_size(x, format="auto")
plotDF <- df %>%
  distinct(TF, `Output type`, file_nrep, .keep_all = TRUE) %>% 
  group_by(`Output type`, file_nrep) %>%
  summarise(
      mean_size = mean(size, na.rm = TRUE),
      n = n()
  )


write_rds(plotDF, str_c(outPrefix, ".n_by_outputType_and_nrep.plotDF.rds"))

p <- plotDF %>%
  ggplot(aes(x = `Output type`, y = n, fill = factor(file_nrep), label = n)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(vjust = 1.2, position = position_dodge(.9)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(str_c(outPrefix, ".n_by_outputType_and_nrep.barplot.pdf"))


p <- plotDF %>%
  ggplot(aes(x = `Output type`, y = mean_size, fill = factor(file_nrep))) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(str_c(outPrefix, ".mean_size_by_outputType_and_nrep.barplot.pdf"))

p <- df %>%
  ggplot(aes(x = `Output type`, y = size, fill = factor(file_nrep))) + 
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(str_c(outPrefix, ".size_by_outputType_and_nrep.boxplot.pdf"))


# Filter for Output tpye tests -------------------------------------------------

fltOuttype <- df %>% 
  arrange(date) %>% 
  filter(TF %in% c("CTCF", "REST", "STAT1", "RAD21")) %>%
  distinct(`Output type`, TF, .keep_all = TRUE) %>%
  mutate(usedTF = TF %in% useTFs) %>%
  arrange(desc(TF), desc(`Output type`)) %>%
  mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`)))

fltOuttype %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.fltOuttype.txt"),
    col_names = FALSE)

fltOuttype %>%
  select(`File accession`:Lab, -rep, filePath) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.fltOuttype.tsv"))

#-------------------------------------------------------------------------------
# Filter for fold change and all TFs -------------------------------------------------
#-------------------------------------------------------------------------------

fcDF <- df %>% 
  # filter for output type "fold change" with 2 replicates or "signal" 
  filter(`Output type` %in% c("fold change over control", "signal")) %>% 
  filter(file_nrep == 2 | `Output type` == "signal") %>%
  # reorder rows to prefere "fold change over control" for each TF in distinct()
  mutate(output_type = factor(`Output type`, c("fold change over control", "signal"))) %>% 
  arrange(TF, output_type) %>% 
  # take only a unique data set per TF
  distinct(TF, .keep_all = TRUE) %>% 
  # annotte with TF selection and file path
  mutate(usedTF = TF %in% useTFs) %>%
  mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`))) %>% 
  select(`File accession`, TF, output_type, usedTF, filePath, file_nrep, Lab, everything())

fcDF %>%
  select(-rep) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.fcDF.tsv"))

fcDF %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.fcDF.txt"),
    col_names = FALSE)

# urls of selected TFs
fcDF %>% 
  filter(TF %in% useTFs) %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.fcDF_selectedTF.txt"),
    col_names = FALSE)

#-------------------------------------------------------------------------------
# HeLa fold-change data
#-------------------------------------------------------------------------------
fcHelaDF <- meta %>%
  filter(`File format` == "bigWig") %>% 
  filter(Assembly == "hg19") %>% 
  filter(`Biosample term name` == "HeLa-S3") %>% 
  # filter for output type "fold change" with 2 replicates or "signal" 
  filter(`Output type` %in% c("fold change over control", "signal")) %>% 
  filter(file_nrep == 2 | `Output type` == "signal") %>%
  # reorder rows to prefere "fold change over control" for each TF in distinct()
  mutate(output_type = factor(`Output type`, c("fold change over control", "signal"))) %>% 
  arrange(TF, output_type) %>% 
  # take only a unique data set per TF
  distinct(TF, .keep_all = TRUE) %>% 
  # annotte with TF selection and file path
  mutate(usedTF = TF %in% useTFs) %>%
  mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`))) %>% 
  select(`File accession`, TF, output_type, usedTF, filePath, file_nrep, Lab, everything())

fcHelaDFselected <- fcHelaDF %>%
  filter(usedTF)
# count(`Biosample term name`)

fcHelaDFselected %>%
  select(-rep) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.fc_HELA_selected.tsv"))

fcHelaDFselected %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.fc_HELA_selected.txt"),
    col_names = FALSE)

