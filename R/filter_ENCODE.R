

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

# "NFKB",
# "NFYB",

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


# Filter for used data sets -----------------------------------------------
# df %>% 
#   filter(`Output type` %in% c("signal")) %>%
#   distinct(TF) %>%
#   summarize(
#     nUsed = sum(TF %in% useTFs),
#     percentUsed = nUsed / length(useTFs) * 100,
#     n = n()
#   )
# 

  # # filter(`Biological replicate(s)` == "1" | is.na(rep)) %>% 
  # # filter(map_lgl(rep, setequal, 1:2) | is.na(rep)) %>% 
  # filter(file_nrep == 2 | `Output type` == "raw signal") %>% 
  # distinct(`Output type`, TF, .keep_all = TRUE) %>%
  # mutate(usedTF = TF %in% useTFs) %>%
  # arrange(desc(usedTF), desc(`Output type`)) %>%
  # mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`)))
  # 

  
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

# Output filtered URL list and metadata table -----------------------------

flt %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.flt.txt"),
    col_names = FALSE)

flt %>%
  select(-rep) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.flt.tsv"))


# Filter for Output tpye tests -------------------------------------------------

fltOuttype <- df %>% 
  # filter(`Output type` %in% c("raw signal", "fold change over control")) %>%
  # filter(`Biological replicate(s)` == "1" | is.na(rep)) %>% 
  # filter(map_lgl(rep, setequal, 1:2) | is.na(rep)) %>% 
  # filter(file_nrep == 2 | `Output type` == "raw signal") %>% 
  # filter(rep == "1" | is.na(rep)) %>%
  # filter(`Biosample term name` == "GM12878") %>%
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

#-------------------------------------------------------------------------------
# get BAM files
#-------------------------------------------------------------------------------

dfBam <- meta %>%
  filter(`File format` == "bam") %>% 
  filter(`Output type` == "alignments") %>% 
  filter(Assembly == "hg19") %>% 
  filter(`Biosample term name` == "GM12878")

# filter for "signal" as output with combined replicate (is NA here) and only GM12878 cell line

fltBam <- dfBam %>% 
  filter(`Biosample term name` == "GM12878") %>%
  distinct(TF, .keep_all = TRUE) %>%
  mutate(usedTF = TF %in% useTFs) %>%
  arrange(desc(usedTF)) %>%
  mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`)))


# Output filtered URL list and metadata table -----------------------------
fltBam %>% 
  select(`File download URL`) %>%
  write_tsv(
    path = file.path("data", "ENCODE", "URLs.fltBam.txt"),
    col_names = FALSE)

fltBam %>%
  select(-rep) %>%
  write_tsv(path = file.path("data", "ENCODE", "metadata.fltBam.tsv"))

#===============================================================================

# sDF <- df %>% 
#   filter(`Output type` %in% c("signal")) %>%
#   # filter(`Biological replicate(s)` == "1, 2") %>%
#   # filter(map_lgl(rep, setequal, 1:2) | is.na(rep)) %>% 
#   # filter(file_nrep == 2 | `Output type` == "raw signal") %>% 
#   # filter(rep == "1" | is.na(rep)) %>%
#   # filter(`Biosample term name` == "GM12878") %>%
#   # filter(TF %in% c("CTCF", "REST", "STAT1")) %>%
#   mutate(usedTF = TF %in% useTFs) %>%
#   arrange(desc(TF), desc(`Output type`), desc(date)) %>%
#   distinct(`Output type`, TF, .keep_all = TRUE) %>%
#   mutate(filePath = file.path("data", "ENCODE", "Experiments", basename(`File download URL`))) %>% 
#   select(1:9, usedTF, filePath, everything())



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


