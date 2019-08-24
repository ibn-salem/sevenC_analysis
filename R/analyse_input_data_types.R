#*******************************************************************************
# Analysis different input data types for loop prediction with sevenC. ------
#*******************************************************************************

library(sevenC)    # devtools::install_github("ibn-salem/sevenC")
library(rtracklayer)  # to import() BED files
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(modelr)       # for tidy modeling
library(precrec)      # for ROC and PRC curves
library(RColorBrewer)   # for nice colors
library(rsample)
library(pryr) # for object_size()
library(feather)      # for efficient storing of data.frames
# library(multidplyr)   # for partition() and collect() to work in parallel
library(furrr)        # for parallelization
library(ROCR)         # for binary clasification metrices

source("R/sevenC.functions.R")

# Set parameter ----------------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(16, parallel::detectCores() - 1)
plan(multicore, workers = N_CORES)

# MIN_MOTIF_SIG <- 5
MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", paste0("v06_input_types.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

DATA_TYPES_META_FILE = "data/DATA_TYPES_metadata_v06.tsv"

# Parse and filter meta data  --------------------------------------------------

meta <- read_tsv(DATA_TYPES_META_FILE, col_types = cols(
  file_accession = col_character(),
  data_type = col_character(),
  TF = col_character(),
  output_type = col_character(),
  path = col_character()
))

meta <- meta %>% 
  mutate(
    name = str_replace_all(str_c(data_type, TF, output_type, sep = "_"), "[ -]", "_"),
    file_exists = file.exists(path),
  ) %>% 
  rename(filePath = path) %>% 
  select(name, file_exists, everything()) %>%
  filter(file_exists)

write_tsv(meta, paste0(outPrefix, ".meta_filtered.tsv"))

# COL_TF <- brewer.pal(8, "Set2")[1:length(unique(meta$TF))]
# names(COL_TF) <- unique(meta$TF)
COL_TF <- c(brewer.pal(8, "Set2")[1:6], brewer.pal(8, "Accent")[2:6])
names(COL_TF) <- unique(meta$TF)
#barplot(1:length(COL_TF), col = COL_TF, names.arg = names(COL_TF))
# pie(rep(1, length(COL_TF)), col = COL_TF, labels = names(COL_TF))

# Select motifs and parse input data -----------------------------------

# read CTCF moitf pairs as candidates
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))
# gi <- read_rds(paste0(
#   str_replace(outPrefix, "v06", "v05"), ".gi.rds"))

# Annotae with coverage and correlation -------------------------------------
# new_idx <- which(str_detect(meta$name, "histone"))

if (!GI_LOCAL ) {
  
  # iterate over all ChIP-seq sata sets
  for (i in seq_len(nrow(meta))) {
  # for (i in new_idx) {
      
    message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
    
    #add coverage and correlation of coverage
    gi <- addCor(
      gi,
      meta$filePath[[i]],
      meta$name[[i]],
    )
  }  
  
  # save file for faster reload
  write_rds(gi, paste0(outPrefix, ".gi.rds"))
  
} else {
  gi <- read_rds(paste0(outPrefix, ".gi.rds"))
}


#*******************************************************************************
# Analyse loopps --------------------------------------------------------
#*******************************************************************************
df <- as_tibble(as.data.frame(mcols(gi))) %>%
  mutate(
    id = 1:nrow(.)
  ) %>% 
  select(id, loop, everything()) 

write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))

#===============================================================================
# Training and prediction in cross-validation
#===============================================================================

# k fold cross validation with 1 repeats
set.seed(3579)

tidyCV <- df %>% 
  vfold_cv(V = K, repeats = 1) %>% 
  tidy()

write_feather(tidyCV, paste0(outPrefix, ".tidyCV.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))


designList <- list(
  "Dist" =  loop ~ dist,
  "Orientation" = loop ~ strandOrientation,
  "Motif" = loop ~ score_min,
  "Dist+Orientation+Motif" = loop ~ dist + strandOrientation + score_min
)

# get design formula for each TF and model
designDF <- tibble(
  meta_name = c(rep(NA, length(designList)), meta$name, meta$name),
  name = c(
    names(designList),
    paste0(meta$name, "_only"),
    meta$name
  ),
  design = c(
    designList,
    map(meta$name, ~as.formula(paste0("loop ~ cor_", .x)) ),
    map(meta$name, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x)) )
  ),
  color = c(
    c("greenyellow", "gold2", "khaki", "brown4"),
    COL_TF[meta$TF],
    COL_TF[meta$TF]
  )
) %>% 
  mutate(name = factor(name, c(paste0(meta$name, "_only"), meta$name, names(designList)))) %>% 
  arrange(name)

write_rds(designDF, paste0(outPrefix, "designDF.rds"))


# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = designDF$name, Fold) %>% 
  # add design formular for each TF
  left_join(designDF, by = "name") %>% 
  mutate(id = parse_integer(str_replace(Fold, "^Fold", "")))



# group data set to clusters
cvDF <- cvDF %>% 
  group_by(name, Fold) %>%
  # fit model on training part
  # fit model and save estimates in tidy format
  mutate(
    tidy_model = future_map2(Fold, design, .f = tidyer_fitter, 
                      tidyCV = tidyCV, data = df),
    pred = pmap(
      list(
        map(Fold, tidy_assessment, data = df, tidyCV = tidyCV),
        design,
        map(tidy_model, "estimate")
      ),
      sevenC:::predLogit
    ),
    label = map(map(Fold, tidy_assessment, data = df, tidyCV = tidyCV), "loop")
  ) %>% 
  ungroup()

write_rds(cvDF, path = paste0(outPrefix, "cvDF.rds"))

