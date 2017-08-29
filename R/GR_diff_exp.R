################################################################################
# Get differentially expressed genes on DEX treatment (from Vocker et al 2016)
################################################################################

require(tidyverse)
require(stringr)
require(DESeq2)

#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------

metadata <- tibble(
  name = c("DEX_3hr_Rep1", "DEX_3hr_Rep2", "EtOh_3hr_Rep1", "EtOH_3hr_Rep2"),
  sample = c("GSM2095194", "GSM2095195", "GSM2095196", "GSM2095197"),
  treatment = factor(c("DEX", "DEX", "EtOH", "EtOH"), c("EtOH", "DEX")),
  replicate = c(1, 2, 1, 2),
  path = c(
    "data/Vockley2016/RNA-seq/GSM2095194_DEX_3hr_Rep1.bam.sorted.counts.txt",
    "data/Vockley2016/RNA-seq/GSM2095195_DEX_3hr_Rep2.bam.sorted.counts.txt",
    "data/Vockley2016/RNA-seq/GSM2095196_EtOH_3hr_Rep1.bam.sorted.counts.txt",
    "data/Vockley2016/RNA-seq/GSM2095197_EtOH_3hr_Rep2.bam.sorted.counts.txt"
  )
)


#-------------------------------------------------------------------------------
# function to preformat each .count file
#-------------------------------------------------------------------------------
parseCountFile <- function(inFile){
  
  # inFile = "data/Vockley2016/RNA-seq/GSM2095194_DEX_3hr_Rep1.bam.sorted.counts.txt"
  
  df_raw <- read_tsv(inFile,  col_names = FALSE, col_type = cols(
    X1 = col_character(),
    X2 = col_integer(),
    X3 = col_integer(),
    X4 = col_integer()
  ))
  
  df <- df_raw %>% 
    mutate(ensg = str_split_fixed(X1, "\\|", 7)[, 2]) %>% 
    mutate(ensg = str_split_fixed(ensg, "\\.", 2)[, 1]) %>% 
    mutate(
      len = X2,
      count = X3
      ) %>% 
    select(ensg, len, count) %>% 
    group_by(ensg) %>% 
    arrange(desc(len)) %>% 
    distinct(ensg, .keep_all = TRUE) %>% 
    select(ensg, count)

  return(df)
}

#-------------------------------------------------------------------------------
# Parse counts and combine to count matrix
#-------------------------------------------------------------------------------
countList <- map(metadata$path, parseCountFile)

countDF <- countList %>%
  Reduce(function(df1, df2) left_join(df1, df2, by = "ensg"), .)

# add column names
names(countDF) <- c("ensg", metadata$name)

# write to output file
write_tsv(countDF, "data/Vockley2016/RNA-seq/countDF.tsv")

#-------------------------------------------------------------------------------
# DESeq2 differentail expression analysis
#-------------------------------------------------------------------------------
# according to http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

dds <- DESeqDataSetFromMatrix(countData = countDF[, 2:ncol(countDF)],
                              colData = metadata,
                              design = ~ treatment,
                              rowData = countDF[, 1])

dds <- DESeq(dds)
rnms <- resultsNames(dds)
res <- results(dds)
# res <- lfcShrink(dds, coef = length(rnms), res = res)

res %>% as.data.frame() %>% 
  as_tibble() %>%
  mutate(ensg = countDF$ensg) %>% 
  dplyr::select(ensg, everything()) %>% 
  write_tsv("data/Vockley2016/RNA-seq/countDF.DESeq2.results.tsv")

