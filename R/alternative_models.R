

require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(colorRamps)  # to pick different colors
require(scales)     # for scales in plotting and percent() formatin function
require(stringr)    # for string function and regular expressions
require(BiocParallel) # for parallelisation
require(randomForest)

# use previously saved gi object?
GI_LOCAL <- FALSE
FACTORBOOK_MOTIFS <- FALSE
CTCF_motif_file <- "data/factorbook/factorbookMotifPos.txt.CTCF.bed"
MIN_MOTIF_SCORE <- 2
OUTPUT_TYPES = FALSE

WINDOW_SIZE <- 1000
BIN_SIZE <- 1
TRUE_LOOPS <- "HIC_ChIAPET"

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path(
  "results", 
  paste0("v01_",
         ifelse(OUTPUT_TYPES, "outTypes", "UCSC"), 
         ifelse(
           FACTORBOOK_MOTIFS, 
           paste0("_factorbook_", MIN_MOTIF_SCORE),
           ""),
         "_w", WINDOW_SIZE,
         "_b", BIN_SIZE
  )
)

ucscMetaFile <- "data/ENCODE_UCSC_FILES.tsv"
outTypeMetaFile <- "data/ENCODE/metadata.fltOuttype.tsv"
# metaFile <- ifelse(OUTPUT_TYPES, outTypeMetaFile, ucscMetaFile)

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

if (OUTPUT_TYPES) {
  # parse metatdata table for input files
  meta <- read_tsv(outTypeMetaFile,
                   col_types = cols(
                     `File accession` = col_character(),
                     TF = col_character(),
                     `Output type` = col_character(),
                     file_nrep = col_character(),
                     exp_nrep = col_integer(),
                     Lab = col_character()
                   )
  )
  
  # add name 
  meta <- meta %>% 
    mutate(name = str_c(TF, "_", str_replace_all(`Output type`, "[ -]", "_"))) %>% 
    select(name, everything())
  
} else {
  # parse ucscMeta file
  meta <- read_tsv(ucscMetaFile,
                   col_types = cols(
                     file = col_character(),
                     TF = col_character()
                   ))
  
  meta <- meta %>% 
    mutate(filePath = file.path("data", "ENCODE", "UCSC", file)) %>% 
    mutate(name = TF)
}


# Load df ----------------------------------------------------------------------
load(paste0(outPrefix, ".df.Rdata"))



useTF <- c(meta$name, "across_TFs")
# useTF <- c("CTCF", "Stat1")

# downsample negative set to the same size as positves
# predDF <- df %>%
#   group_by(loop) %>%
#   sample_n(nLoop) %>%
#   ungroup()

# predDF <- df 
# %>% 
#   filter(strandOrientation == "convergent")



# dived in training and test set by taking 9/10 for training and 1/10 for testing
testCases <- sample.int(nrow(df), size = round(nrow(df)/10))
train <- df[-testCases,]
test <- df[testCases,] %>% 
  mutate(strandOrientation = parse_factor(strandOrientation, levels = NULL))

nLoop <- sum(train$loop == "Loop")
# downsample negative set to the same size as positves
trainDS <- train %>% 
  group_by(loop) %>%
  sample_n(nLoop) %>%
  ungroup()

# get predictions from single TF with dist

# for (name in ucscMeta$name) {
# for (name in useTF) {
# predList <- bplapply(useTF, function(name) {
predList <- lapply(useTF, function(name) {
  
  message("INFO: Fit model for TF: ", name)
  
  # design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", "cor_", str_replace(name, "-", "_")))
  # design <- as.formula(paste0("loop ~ ", " dist + cor_", name))
  
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + cor_", name))
  
  # model <- glm(
  #   design, 
  #   family = binomial(link = 'logit'), 
  #   data = train)
  
  trainSub <- train %>% 
    select(
      loop, 
      dist,
      strandOrientation,
      one_of(str_c("cor_", name))
    ) %>% 
    mutate(strandOrientation = parse_factor(strandOrientation, levels = NULL)) %>% 
    na.omit %>% 
    sample_n(10^5)

  
  fit <- randomForest(
    formula = design, 
    data = trainSub)
  
  predictions <- predict(fit, test)
  
  # test[, str_c("pred_", name)] <- predict(model, newdata = test, type = 'response')
  return( predict(model, newdata = test, type = 'response') )
})
names(predList) <- str_c("pred_", useTF)

predDSList <- lapply(useTF, function(name) {
  
  message("INFO: Fit model for TF: ", name)
  design <- as.formula(paste0("loop ~ ", " dist + strandOrientation + cor_", name))
  
  modelDS <- glm(
    design, 
    family = binomial(link = 'logit'), 
    data = trainDS)
  
  message("INFO: Predict loops with model: ", name)
  
  # test[, str_c("pred_", name)] <- predict(model, newdata = test, type = 'response')
  return( predict(modelDS, newdata = test, type = 'response') )
})
names(predDSList) <- str_c("pred_", useTF, "_DS")

predDF <- as_tibble(c(predList, predDSList))

test <- bind_cols(test, predDF)

# add modles of all TF and only dist

designAll <- as.formula(paste0("loop ~ dist + strandOrientation + ", paste(paste0("cor_", useTF), collapse = " + ")))
modelAll <- glm(
  designAll,
  family = binomial(link = 'logit'),
  data = train)

designDist <- as.formula(paste0("loop ~ dist + strandOrientation"))
# designDist <- as.formula(paste0("loop ~ dist"))
modelDist <- glm(
  designDist, 
  family = binomial(link = 'logit'), 
  data = train)

test[, "pred_all_TF"] <- predict(modelAll, newdata = test, type = 'response')
test[, "pred_dist"] <- predict(modelDist, newdata = test, type = 'response')

# test <- test %>% 
#   mutate(
#     pred_all <- predict(modelAll, newdata = ., type = 'response'),
#     pred_dist <- predict(modelDist, newdata = ., type = 'response')
#   )

save(train, test, file = paste0(outPrefix, ".train_test.Rdata"))

