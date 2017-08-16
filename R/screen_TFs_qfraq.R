################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(modelr)       # for tidy modeling
require(forcats)      # for factors
require(precrec)      # for ROC and PRC curves
require(RColorBrewer)   # for nice colors
require(rtracklayer)  # to import() BED files
require(colorRamps)  # to pick different colors
require(scales)     # for scales in plotting and percent() formatin function
require(stringr)    # for string function and regular expressions
require(BiocParallel) # for parallelisation
require(biobroom)     # to make BioC classes tidy
require(grid)         # for textGrob
require(magrittr)     # for "tee" pipe %T>% 
require(rsample)
# require(recipes)
require(pryr) # for object_size()

#-----------------------------------------------------------------------
# Options for parallel computation
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed = 34312)
# set options
register(multicorParam)
# bpparam() # to print current options


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
MIN_MOTIF_SIG <- 6


WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation

TRUE_LOOPS <- "HIC_ChIAPET"
# TRUE_LOOPS <- "HIC_ChIAPET_CaptureC"

# COL_TF = c(colorRampPalette(brewer.pal(8, "Set1"))(9), "#80da3a")
# COL_TF = c(colorRampPalette(brewer.pal(12, "Set3"))(10), "gray70", "gray50", "gray30")
# pie(rep(1, length(COL_TF)), col=COL_TF, labels=COL_TF, main=length(COL_TF))

# COL_TF = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# str_view(grDevices::colors(), "[^\d]$")
# str_view(grDevices::colors(), "^((?!gr(a|e)y).)*[^\\d]$")
COL_TF <- grDevices::colors()[str_detect(grDevices::colors(), "^((?!(gr(a|e)y|white)).)*[^\\d]$")]
COL_TF <- COL_TF[!COL_TF %in% c("white", "black", "aliceblue", "azure", "beige", "bisque", "cornsilk", "cyan", "darkorchid", "coral", "darkmagenta")]
# pie(rep(1, 50), col=COL_TF, labels = COL_TF)
# pie(rep(1, length(COL_TF)), col=COL_TF)
# COL_TF = colorRamps::primary.colors(12)

#barplot(1:10, col=COL_TF)
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

# outPrefix <- file.path("results", "v01_UCSC_motifs")
outPrefix <- file.path(
  "results", 
  paste0("v03_screen_TF_qfraq.",
         paste0("motifSig", MIN_MOTIF_SIG),
         "_w", WINDOW_SIZE,
         "_b", BIN_SIZE
  )
)

dir.create(dirname(outPrefix), showWarnings = FALSE)

# True loops in GM12878 from Rao et al:
LoopRao2014_GM12878_File <- 
  "data/Rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt"

# ChIA-PET loops in GM12878 from Tang et al 2015:
LoopTang2015_GM12878_Files <- c(
  "data/Tang2015/GSM1872886_GM12878_CTCF_PET_clusters.txt",
  "data/Tang2015/GSM1872887_GM12878_RNAPII_PET_clusters.txt")

CaptureHiC_Files <- c(
  "data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt",
  "data/Mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt"
)

metaFile <- "data/ENCODE/metadata.fltBam.tsv"

#-------------------------------------------------------------------------------
# Parse and filter input ChiP-seq data  -----------------------------------
#-------------------------------------------------------------------------------

# parse ucscMeta file
meta <- read_tsv(metaFile,
                 col_types = 
                   cols(
                     `File accession` = col_character(),
                     TF = col_character(),
                     `Output type` = col_character(),
                     file_nrep = col_character(),
                     exp_nrep = col_integer(),
                     Lab = col_character(),
                     filePath = col_character()
                   )
)

meta <- meta %>% 
  mutate(name = TF) %>%
  mutate(filePath = str_c(filePath, "-qfrags_allChr_chip.bed.sorted.bedGraph.bw")) %>% 
  # filter(file.exists(filePath)) %>% 
  select(TF, name, filePath, everything())

# meta <- meta[1:3, ]

# Select motifs and parse input data -----------------------------------
if (!GI_LOCAL) {
  
  
  ancGR <- chromloop::motif.hg19.CTCF
  
  # filter for p-valu <= MIN_MOTIF_SIG
  ancGR <- ancGR[ancGR$sig >= MIN_MOTIF_SIG]
  
  seqInfo <- seqinfo(chromloop::motif.hg19.CTCF)
  
  # get all pairs within 1M distance
  gi <- chromloop::getCisPairs(ancGR, maxDist = 10^6)
  
  # add strand combinations
  gi <- chromloop::addStrandCombination(gi)
  
  # add motif score
  gi <- chromloop::addMotifScore(gi, colname = "sig")
  
  # parse loops
  trueLoopsRao <- chromloop::parseLoopsRao(
    LoopRao2014_GM12878_File, seqinfo = seqInfo)
  trueLoopsTang2015 <- do.call(
    "c",
    lapply(LoopTang2015_GM12878_Files, 
           chromloop::parseLoopsTang2015, 
           seqinfo = seqInfo))
  
  trueCaptureHiC <- do.call(
    "c",
    lapply(CaptureHiC_Files,
           chromloop::parseCaptureHiC, seqinfo = seqInfo)
  )
  trueCaptureHiC <- trueCaptureHiC[trueCaptureHiC$log_observed_expected >= 10]
  # trueCaptureHiC <- chromloop::parseCaptureHiC(
  #   CaptureHiC_Files[1], seqinfo = seqInfo)
  
  # ol <- IRanges::overlapsAny(gi, trueLoops)
  # gi$Loop_Rao_GM12878 <- factor(ol, c(FALSE, TRUE), c("No loop", "Loop"))
  
  gi <- addInteractionSupport(gi, trueLoopsRao, "Loop_Rao_GM12878")
  gi <- addInteractionSupport(gi, trueLoopsTang2015, "Loop_Tang2015_GM12878")
  gi <- addInteractionSupport(gi, trueCaptureHiC, "Loop_Mifsud2015_GM12878")
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.tmp.Rdata"))
  
}else{
  load(paste0(outPrefix, ".gi.tmp.Rdata"))  
}

# Annotae with coverage and correlation -------------------------------------

if (!GI_LOCAL ) {
  
  # iterate over all ChIP-seq sata sets
  for (i in seq_len(nrow(meta))) {
    
    # check if sample is not present yet.
    if ( !paste0("cor_", meta$name[i]) %in% names(mcols(gi))){
      
      message("INFO: --> Working on sample: ", meta$name[i], ", ", i, " of ", nrow(meta), " <--")
      
      # add coverage  
      regions(gi) <- chromloop::addCovToGR(
        regions(gi), 
        meta$filePath[i], 
        window = WINDOW_SIZE,
        bin_size = BIN_SIZE,
        colname = paste0("cov_", meta$name[i])
      )
      
      # add correlations
      gi <- chromloop::applyToCloseGI(
        gi, 
        datcol = paste0("cov_", meta$name[i]),
        fun = cor, 
        colname = paste0("cor_", meta$name[i])
      )  
    }
    
  }
  
  # Annotae with correlation across TFs -------------------------------------
  
  # function to sum up all vectors in list
  sumList <- function(l){
    map_dbl(l, sum)
  }
  
  covCols <- paste0("cov_", meta$name)
  covDF <- as_tibble(as.data.frame(mcols(regions(gi))[,covCols])) %>% 
    dplyr::mutate_all(.funs = sumList)
  
  mcols(regions(gi))[, "cov_sum"] <- NumericList(as_tibble(t(covDF)))
  
  gi <- chromloop::applyToCloseGI(
    gi, 
    datcol = "cov_sum",
    fun = cor, 
    colname = "cor_across_TFs"
  )
  
  # save file for faster reload
  save(gi, file = paste0(outPrefix, ".gi.Rdata"))
  
} else {
  load(paste0(outPrefix, ".gi.Rdata"))  
}


#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------

df <- as_tibble(as.data.frame(mcols(gi))) %>%
  mutate(
    id = 1:nrow(.),
    HIC_ChIAPET_CaptureC = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop" | Loop_Mifsud2015_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop")),
    HIC_ChIAPET = factor(
      Loop_Tang2015_GM12878 == "Loop" | Loop_Rao_GM12878 == "Loop",
      c(FALSE, TRUE),
      c("No loop", "Loop"))
  ) 


if (TRUE_LOOPS == "HIC_ChIAPET_CaptureC") {
  df <- df %>% mutate(loop =  HIC_ChIAPET_CaptureC)
  outPrefix <- str_c(outPrefix, "_HIC_ChIAPET_CaptureC")
} else {
  df <- df %>% mutate(loop =  HIC_ChIAPET)
}

df <- df %>% 
  select(id, loop, everything()) 

save(df, file = paste0(outPrefix, ".df.Rdata"))  
# load(paste0(outPrefix, ".df.Rdata"))

# remove all but 3 TF columns
rmNames <- paste0("cor_", meta$name[4:nrow(meta)])
df <- df %>%
  select(-match(rmNames, names(.)))

# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = name, value = cor) %>% 
  mutate(name = str_sub(name, 5))


# ------------------------------------------------------------------------------
# Compare Correlation with Boxplot
# ------------------------------------------------------------------------------

tidySubDF <- tidyDF %>% 
  sample_n(min(nrow(tidyDF), 10^7))

p <- ggplot(tidySubDF, aes(x = name, y = cor, col = loop)) +
  geom_boxplot() +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    # text = element_text(size=20), 
    # axis.text.x=element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Correlation of ChIP-seq signal", x="")
ggsave(p, file = paste0(outPrefix, ".cor.by_TF_and_loop.boxplot.pdf"), w = 28, h = 7)

# ------------------------------------------------------------------------------
# Predict loops 
# ------------------------------------------------------------------------------

useTF <- meta$name

#filter for a subset of TFS
useTF <- useTF[1:3]


# ------------------------------------------------------------------------------
# cross validation in a tidy approach 
# ------------------------------------------------------------------------------

# k fold cross validation with 1 repeats
set.seed(3579)

dfCV <- df %>% 
  vfold_cv(V = K, repeats = 1)

# get design formula for each TF
designDF <- tibble(
  name = useTF,
  design = map(useTF, ~as.formula(paste0("loop ~ dist + strandOrientation + score_min + cor_", .x)) )
)

#' Get prediciton result on test subset
#' Source: https://github.com/topepo/rsample/blob/master/vignettes/Working_with_rsets.Rmd
#' 
#' @param splits an individual \code{rsplit} object.
#' 
holdout_results <- function(splits, formula, dependant_var = "loop", ...) {
  
  # Fit the model to the training data set
  mod <- glm(formula = formula, 
             data = analysis(splits), 
             family = binomial(),
             model = FALSE,  
             x = FALSE,
             ...)
  
  # `augment` will save the predictions with the holdout data set
  res <- broom::augment(
      mod, 
      newdata = assessment(splits), 
      type.predict = "response"
    ) %>% 
    as_tibble() %>% 
    select(label = one_of(dependant_var), pred = ".fitted")
  

  # Return the assessment data set with the additional columns
  return(res)
}

# expand data.frame to have all combinations of model and split
dfCV <- dfCV %>% 
  expand(name = useTF, id) %>% 
  left_join(dfCV, by = "id") %>% 
  # add design formular for each TF
  left_join(designDF, by = "name")


# fit model and add predictions for each combination of model and split
dfCV <- dfCV %>% 
  mutate(
    result = map2(splits, design, holdout_results)
  )

# extract only the needed columns
evalDF <- dfCV %>% 
  mutate(
    pred = map(result, "pred"),
    label = map(result, "label")
  ) %>% 
  mutate(
    fold = parse_integer(str_replace(id, "Fold", ""))
  ) %>% 
  select(name, fold, pred, label)

# save evalDF
save(evalDF, file = paste0(outPrefix, ".evalDF.Rdata"))  
# load(paste0(outPrefix, ".evalDF.Rdata"))


# get number and percent of positives
posDF <- evalDF %>% 
  mutate(
    n_pos = map_dbl(label, function(l) sum(l == "Loop", na.rm = TRUE)),
    percent_pos = n_pos / map_dbl(label, length) * 100
  )


# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = evalDF$pred,
  labels = evalDF$label,
  modnames = evalDF$name,
  dsids = evalDF$fold,
  posclass = levels(evalDF$label[[1]])[2],
  x_bins = 100)


# get data.frame with auc values
aucDF <- as_tibble(auc(curves))

# get ranked modle names
ranked_models <- aucDF %>% 
  filter(curvetypes == "PRC") %>% 
  group_by(modnames) %>% 
  summarize(
    auc_mean = mean(aucs, na.rm = TRUE)
  ) %>% 
  arrange(auc_mean) %>% 
  select(modnames) %>% 
  unlist()

# order aucDF by ranks
aucDF <- aucDF %>% 
  mutate(
    modnames = factor(modnames, ranked_models)
  ) %>% 
  arrange(desc(modnames))

# get data from fro ggplot
# curveDF <- precrec::fortify(curves)

# build color vector with TFs as names
# TF_models <- ranked_models
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- ranked_models

###################
# TODO Continue here!
###################

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)

# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curves, "ROC") + 
  # coord_cartesian(expand = FALSE) + 
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFroc$aucs,3)
                     ),
                     guide = guide_legend(
                       override.aes = list(size = 2),
                       reverse = TRUE)
  ) + 
  theme(legend.position = c(.75,.4))

# g
ggsave(g, file= paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(curves, "PRC", size = 4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFprc$modnames, 
                       ": AUC=", 
                       signif(aucDFprc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position = "none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)
#===============================================================================
#===============================================================================


# ------------------------------------------------------------------------------
# New tidy approach :)
# ------------------------------------------------------------------------------

design <- loop ~ dist + strandOrientation + score_min + cor

grp <- tidyDF %>% 
  select(-starts_with("Loop_"), -score_1, -score_2) %>% # remove unused columns
  filter(name %in% useTF) %>% # filter for only TFs in useTF
  group_by(name)              # group by TF

# Make k-fold cros-validation by splitting data for each gorup (TF)
cv <- grp %>% 
  do(crossv_kfold(data = ., k = K)) 


#' fitts a logit model to training data and returns prediction response vector
#' on test data
#' 
#' @param training_data data set with training data.
#' @param design a formular object defining the prediction design and variables.
fitter <- function(training_data, design) {
  glm(
    design, 
    family = binomial(link = 'logit'),
    data = training_data,
    model = FALSE, 
    x = FALSE
  ) 
}


#' fitts a logit model to training data and returns prediction response vector
#' on test data
#' 
#' @param training_data data set with training data.
#' @param test_data data set with training data.
#' @param design a formular object defining the prediction design and variables.
fitter_augment <- function(training_data, test_data, design) {
  glm(
    design, 
    family = binomial(link = 'logit'),
    data = training_data,
    model = FALSE, 
    x = FALSE
  ) %>% 
  broom::augment(newdata = training_data, type.predict = "response") %>% 
  as_tibble()
}


#' fitts a logit model to training data and returns model fit
#' 
#' @param training_data data set with training data.
#' @param design a formular object defining the prediction design and variables.
fitter <- function(training_data, design) {
  glm(
    design, 
    family = binomial(link = 'logit'),
    data = training_data,
    model = FALSE, 
    x = FALSE
  ) 
}

# %>% 
#   predict(object = ., newdata = testing_data, type = "response")

require(biglm)
bm <- biglm::bigglm(design, data = cv$train[[1]], family = binomial(link = 'logit'))

# apply traingin and prediction to data sets
cvDF <- cv %>% 
  # mutate(
  #   model = map(train, fitter, design = design)
  # ) %>% 
  mutate(
    aug = map2(train, test, fitter_augment, design = design)
  )
# ,
#     pred = map2(model, test, prediction, type = "response"),
#     label = map(test, function(df) as_tibble(df)[["loop"]])
#   ) %>% 
#   select(-train, -test) %>% 
#   mutate(fold = parse_integer(.id))

# pryr::object_size(tidyDF)
save(cvDF, file = paste0(outPrefix, ".cvDF.Rdata"))  
# load(paste0(outPrefix, ".cvDF.Rdata"))

#-------------------------------------------------------------------------------
# Analyse Model 
#-------------------------------------------------------------------------------

cvDF <- cvDF %>% 
  mutate(
    n_pos = map_dbl(label, function(l) sum(l == "Loop", na.rm = TRUE))
  ) %>% 
  mutate(tidy_model = map(model, broom::tidy)) %>% 
  mutate(glance = map(model, broom::glance))


# add model quality and mdoel as tidy DF
byTFfold <- byTFfold %>% 
  mutate(tidy_model = map(model, broom::tidy)) %>% 
  mutate(glance = map(model, broom::glance))

save(byTFfold, file = paste0(outPrefix, ".byTFfold.Rdata"))

# Analyse Model parameters ----------------------------------------------------

# get DF with model quality
modelQualDF <- byTFfold %>% 
  unnest(glance, .drop = TRUE)
write_tsv(modelQualDF, paste0(outPrefix, ".modelQualDF.tsv"))

# get tidy model DF
modelDF <- byTFfold %>% 
  unnest(tidy_model, .drop = TRUE)

write_tsv(modelDF, paste0(outPrefix, ".modelDF.tsv"))

# add meta data
modelDF <- modelDF %>% 
  left_join(meta, by = "name")

paramByTF <- modelDF %>% 
  group_by(name, term) %>% 
  summarize(
    n = n(),
    estimate_mean = mean(estimate),
    estimate_sd = sd(estimate)
  )

p <- ggplot(paramByTF, aes(x = name, y = estimate_mean, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(
    ymin = estimate_mean - estimate_sd, 
    ymax = estimate_mean + estimate_sd), width = 0.25) +
  # geom_text(aes(label = round(estimate, 2)), vjust = "inward") + 
  # facet_grid(param ~ ., scales = "free_y") + 
  geom_text(aes(label = round(estimate_mean, 2)), hjust = "inward") + 
  facet_grid(. ~ term , scales = "free_x") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, file = paste0(outPrefix, ".paramter.barplot.pdf"), w = 6, h = 12)




#====================
#====================
#====================



# tidy prediction --------------------------------------------------
byTF <- tidyDF %>%
  filter(name != "across_TFs") %>% 
  filter(name %in% useTF) %>% 
  mutate(fold = rep(fold, length(useTF))) %>% 
  group_by(name) %>% 
  nest()


save(byTF, file = paste0(outPrefix, ".byTF_TMP.Rdata"))  
# load(paste0(outPrefix, ".byTF_TMP.Rdata"))


#' Add prediciton to a data set using a glm model and design
#' 
add_pred <- function(df, design = loop ~ dist + strandOrientation + score_min + cor)){
  
  folds <- unique(df$fold)
  
  cvDF <- tibble(
    folds = folds,
    df = list(df)
  )
  
  cvDF <- cvDF %>% 
    mutate(
      train = map2(df, folds, function(d, k) filter(d, fold != k)),
      test = map2(df, folds, function(d, k) filter(d, fold == k))
    ) %>% 
    mutate(
      pred = map2(train, test, function(trainDF, testDF){
        predict(
          object = glm(
            design, 
            family = binomial(link = 'logit'), 
            data = trainDF),
          newdata = testDF,
          type = "response")
      })
    )
  
  model <- map(folds, function(k) {
    glm(
      design, 
      family = binomial(link = 'logit'), 
      data = subset(df, fold != k)
    )
  })
  
  cvDF <- cvDF %>% 
    mutate(
      id = map(test, "id"),
      idpred = map2(id, pred, tibble)
    ) %>% 
    select(idpred)
    
  
  # add prediction to df
  cvDF <- cvDF %>% 
    mutate(
     
    )    
}

singleTF_model <- function(df, fold_idx = 1) {
  
  glm(
    loop ~ dist + strandOrientation + score_min + cor, 
    family = binomial(link = 'logit'), 
    data = subset(df, fold != fold_idx)
  )  
}

add_predict_fold <- function(df, model, fold_idx = 1) {
  df[df$fold == fold_idx, "pred"] <- predict(
    model, 
    newdata = subset(df, fold == fold_idx), 
    type = "response")
  return(df)
}

add_only_predict_fold <- function(df, fold_idx = 1) {

  model <- glm(
    loop ~ dist + strandOrientation + score_min + cor, 
    family = binomial(link = 'logit'), 
    data = subset(df, fold != fold_idx)
  )  
  
  df[df$fold == fold_idx, "pred"] <- predict(
    model, 
    newdata = subset(df, fold == fold_idx), 
    type = "response")
  return(df)
}


# m <- singleTF_model(byTF$data[[1]], 2)
# models <- map(byTF$data, singleTF_model)

for (k in 1:K) {
  
  message("INFO: k = ", k)
  
  # byTF <- byTF %>%
  #   mutate(
  #     !!paste0("model_", k) := map(byTF$data, singleTF_model, fold_idx = k)
  #   )
  
  byTF <- byTF %>% 
    mutate(
      data = map(data, .f = add_only_predict_fold, fold_idx = k)
    )
}

# save(byTF, file = paste0(outPrefix, ".byTF_TESTING.Rdata"))
save(byTF, file = paste0(outPrefix, ".byTF.Rdata"))
# load(paste0(outPrefix, ".byTF.Rdata"))
#-------------------------------------------------------------------------------
# unpack to by TF and fold
#-------------------------------------------------------------------------------

# gather models by fold
byTFfold <- byTF %>%
  gather(starts_with("model_"), key = "fold", value = "model") %>% 
  mutate(
    fold = parse_integer(str_replace(fold, "model_", "")), 
    data = map2(data, fold, function(df, k) filter(df, fold == k))
    )

# add model quality and mdoel as tidy DF
byTFfold <- byTFfold %>% 
  mutate(tidy_model = map(model, broom::tidy)) %>% 
  mutate(glance = map(model, broom::glance))
  
save(byTFfold, file = paste0(outPrefix, ".byTFfold.Rdata"))

# Analyse Model parameters ----------------------------------------------------

# get DF with model quality
modelQualDF <- byTFfold %>% 
  unnest(glance, .drop = TRUE)
write_tsv(modelQualDF, paste0(outPrefix, ".modelQualDF.tsv"))

# get tidy model DF
modelDF <- byTFfold %>% 
  unnest(tidy_model, .drop = TRUE)

write_tsv(modelDF, paste0(outPrefix, ".modelDF.tsv"))

# add meta data
modelDF <- modelDF %>% 
  left_join(meta, by = "name")

paramByTF <- modelDF %>% 
  group_by(name, term) %>% 
  summarize(
    n = n(),
    estimate_mean = mean(estimate),
    estimate_sd = sd(estimate)
  )

p <- ggplot(paramByTF, aes(x = name, y = estimate_mean, fill = name)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(
    ymin = estimate_mean - estimate_sd, 
    ymax = estimate_mean + estimate_sd), width = 0.25) +
  # geom_text(aes(label = round(estimate, 2)), vjust = "inward") + 
  # facet_grid(param ~ ., scales = "free_y") + 
  geom_text(aes(label = round(estimate_mean, 2)), hjust = "inward") + 
  facet_grid(. ~ term , scales = "free_x") + 
  coord_flip() +
  labs(y = "Parameter estimate", x = "Model") + 
  theme_bw() + scale_fill_manual(values = COL_TF) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, file = paste0(outPrefix, ".paramter.barplot.pdf"), w = 6, h = 12)

#-------------------------------------------------------------------------------
# Analyse performace 
#-------------------------------------------------------------------------------

evalDF <- byTFfold %>% 
  mutate(
    pred =  map(data, function(df) df[["pred"]]),
    label =  map(data, function(df) df[["loop"]])
  ) %>% 
  select(name, fold, pred, label)

posDF <- evalDF %>% 
  mutate(
    n_pos = map_dbl(label, function(l) sum(l == "Loop", na.rm = TRUE))
  )

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
      scores = evalDF$pred,
      label = evalDF$label,
      modnames = evalDF$name,
      dsids = evalDF$fold,
      calc_avg = TRUE)

# get data.frame with auc values
aucDF <- as_tibble(auc(curves))

# get ranked modle names
ranked_models <- aucDF %>% 
  filter(curvetypes == "PRC") %>% 
  group_by(modnames) %>% 
  summarize(
    auc_mean = mean(aucs, na.rm = TRUE)
  ) %>% 
  arrange(auc_mean) %>% 
  select(modnames) %>% 
  unlist()

# order aucDF by ranks
aucDF <- aucDF %>% 
  mutate(
    modelnames = factor(modelnames, ranked_models)
  ) %>% 
  arrange(desc(modelnames))

# get data from fro ggplot
curveDF <- precrec::fortify(curves)

# build color vector with TFs as names
# TF_models <- ranked_models
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- ranked_models

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)

# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curves, "ROC") + 
  # coord_cartesian(expand = FALSE) + 
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFroc$aucs,3)
                      ),
                     guide = guide_legend(
                      override.aes = list(size = 2),
                      reverse = TRUE)
                     ) + 
  theme(legend.position = c(.75,.4))

# g
ggsave(g, file= paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(curves, "PRC", size = 4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFprc$modnames, 
                       ": AUC=", 
                       signif(aucDFprc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position = "none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)



autoplot(curves, show_cb = TRUE)
#==============================================
# OLD CODE BELLOW
#==============================================

# model_names <- c(useTF, paste0(useTF, "_DS"), "dist", "all_TF")
model_names <- useTF

# plot ROC and PRC
modelData <-  mmdata(
  # scores = as.list(select(test, one_of(str_c("pred_", model_names)))),
  scores = map(byTF$data, function(subDF) filter(subDF, fold == 1)[["pred"]]),
  labels = filter(df, fold == 1)[["loop"]],
  modnames = model_names,
  posclass = "Loop"
)

# get rank of models
ranked_models <- auc( evalmod(modelData) ) %>% 
  filter(curvetypes == "PRC") %>% 
  arrange(aucs) %>% 
  select(modnames) %>% 
  unlist()

# build color vector with TFs as names
# TF_models <- ranked_models
COL_TF <- COL_TF[seq(1, length(ranked_models))]
names(COL_TF) <- ranked_models


# build model again with ordered modelnames
modelData <-  mmdata(
  scores = map(byTF$data, function(subDF) filter(subDF, fold == 1)[["pred"]]),
  labels = filter(df, fold == 1)[["loop"]],
  modnames = ranked_models,
  posclass = "Loop"
)

# caluclate ROC, PRC, and basic performance scores
curves <- evalmod(modelData)
# scores <- evalmod(modelData, mode = "basic")

# get AUC of ROC and PRC
aucDF <- as_tibble(auc(curves)) %>%
  left_join(meta, by = c("modnames" = "name")) %>% 
  mutate(modnames = fct_reorder(modnames, aucs, fun = min)) %T>% 
  write_tsv(paste0(outPrefix, ".aucDF.tsv"))

COL_TF <- COL_TF[1:length(model_names)]
names(COL_TF) <- levels(aucDF$modnames)


# %>%
#   mutate(modnames = factor(modnames, ranked_models))

# barplot of AUCs of ROC and PRC

p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)

p <- ggplot(aucDF, aes(x = TF, y = aucs, fill = modnames)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(aucs, 2)), size = 3, hjust = 1, angle = 90) +
  facet_grid(curvetypes ~ outType, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = COL_TF) +
  labs(x = "Models", y = "AUC")
# p
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.by_TF.barplot.pdf"), w = 14, h = 7)


# get ROC plots
aucDFroc <- aucDF %>% 
  filter(curvetypes == "ROC")

g <- autoplot(curves, "ROC") + 
  # coord_cartesian(expand = FALSE) + 
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFroc$modnames, 
                       ": AUC=", 
                       signif(aucDFroc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) + 
  theme(legend.position=c(.75,.4))
# g
ggsave(g, file= paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)

# get PRC plots
aucDFprc <- aucDF %>% 
  filter(curvetypes == "PRC")

g <- autoplot(curves, "PRC", size = 4) +
  scale_color_manual(values = COL_TF,
                     labels = paste0(
                       aucDFprc$modnames, 
                       ": AUC=", 
                       signif(aucDFprc$aucs,3)),
                     guide = guide_legend(override.aes = list(size = 2),
                                          reverse = TRUE)) +
  # theme(legend.position=c(.75,.6))
  theme(legend.position = "none")
# g
ggsave(g, file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)

