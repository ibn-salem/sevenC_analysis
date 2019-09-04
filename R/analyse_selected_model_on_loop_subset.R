################################################################################
# Analysis of predictd chromatin looping interactions using the sevenC tool
################################################################################


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
library(ROCR)         # for binary clasification metrices
library(furrr)        # for parallelization

source("R/sevenC.functions.R")

# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(10, parallel::detectCores() - 1)

MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

# use mutliple cores
plan(multicore, workers = N_CORES)

CTCF_peaks_bed <- "data/ENCODE/Peaks/ENCFF710VEH.bed"

outPrefixSelectedModel <- file.path("results", paste0("v05_selected_models.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

outPrefix <- file.path("results", paste0("v05_selected_models_on_loop_subset.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

screenPrefix <- file.path("results", paste0("v05_screen_TF_lfc.", 
                                            paste0("motifPval", MOTIF_PVAL), 
                                            "_w", WINDOW_SIZE, 
                                            "_b", BIN_SIZE))

# metadata file
metaFile <- "data/ENCODE/metadata.fcDF.tsv"

SELECTED_TF <- c(
  "RAD21",
  "CTCF",
  "ZNF143",
  "STAT1",
  "EP300",
  "POLR2A"
)

COL_SELECTED_TF = brewer.pal(length(SELECTED_TF), "Set1")
COL_SELECTED_TF_1 = brewer.pal(12, "Paired")[c(1, 3, 5, 7, 9, 11)]
COL_SELECTED_TF_2 = brewer.pal(12, "Paired")[c(2, 4, 6, 8, 10, 12)]
names(COL_SELECTED_TF_1) <- SELECTED_TF
names(COL_SELECTED_TF_2) <- SELECTED_TF

LOOP_LEVEL <- c("No loop", "Loop")
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

# reformat metadata
meta <- meta %>% 
  filter(TF %in% SELECTED_TF) %>% 
  # mutate(name = paste0(TF, "_lfc")) %>%
  mutate(name = TF) %>% 
  arrange(factor(name, SELECTED_TF)) %>% 
  select(TF, name, filePath, everything())


write_tsv(meta, paste0(outPrefix, ".meta.tsv"))
# meta <- read_tsv(paste0(outPrefix, ".meta.tsv"))

# Select motifs and parse input data -----------------------------------

# read CTCF moitf pairs with ChIP-seq data

# gi <- read_rds(paste0(outPrefixSelectedModel, ".gi.rds"))

# # read CTCF moitf pairs as candidates
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))


#-------------------------------------------------------------------------------
# Analyse loopps --------------------------------------------------------
#-------------------------------------------------------------------------------
df <- read_feather(paste0(outPrefixSelectedModel, ".df.feather"))

# add union and intersection
df <- df %>% 
  mutate(
    Loop_Union = loop,
    Loop_Intersection = factor(
      Loop_Rao_GM12878 == "Loop" & 
      Loop_Tang2015_GM12878_CTCF == "Loop" &
      Loop_Tang2015_GM12878_RNAPII == "Loop",
      c(FALSE, TRUE), 
      c("No loop", "Loop")
    )
  )

#===============================================================================
# Training and prediction in cross-validation
#===============================================================================

# k fold cross validation with 1 repeats
set.seed(3579)

tidyCV <- read_feather(paste0(outPrefixSelectedModel, ".tidyCV.feather"))

dependent_variable <-  c("Loop_Rao_GM12878", "Loop_Tang2015_GM12878_CTCF", 
                         "Loop_Tang2015_GM12878_RNAPII", "Loop_Union", "Loop_Intersection")

designDF <- tibble(
  dep_var = dependent_variable,
  design = map(dep_var, function(dep_var) 
    c(
      list(
        "Dist" =  as.formula(str_c(dep_var, " ~ dist")),
        "Orientation" =  as.formula(str_c(dep_var, " ~ strandOrientation")),
        "Motif" =  as.formula(str_c(dep_var, " ~ score_min")),
        "Dist+Orientation+Motif" =  as.formula(str_c(dep_var, " ~ dist + strandOrientation + score_min"))
      ),
      map(meta$name, ~as.formula(str_c(dep_var, " ~ dist + strandOrientation + score_min + cor_", .x)) )
    )
  ),
  name = list(c("Dist", "Orientation", "Motif", "Dist+Orientation+Motif", meta$name)),
  color = list(c("greenyellow", "gold2", "khaki", "brown4", COL_SELECTED_TF_2))
) %>% 
  unnest(design, name, color) %>% 
  mutate(
    dep_var = factor(dep_var, unique(dep_var)),
    name = factor(name, unique(name))
    ) %>% 
  arrange(dep_var, name)

write_rds(designDF, paste0(outPrefix, "designDF.rds"))

# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = designDF$name, Fold) %>% 
  # add design formular for each TF
  left_join(designDF, by = "name") %>% 
  mutate(id = parse_integer(str_replace(Fold, "^Fold", "")))

if (!GI_LOCAL){
  # fit model on training part
  cvDF <- cvDF %>% 
    group_by(name, Fold) %>% 
    # fit model and save estimates in tidy format
    mutate(
      tidy_model = future_map2(Fold, design, .f = tidyer_fitter, 
                        tidyCV = tidyCV, data = df),
      pred = future_pmap(
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

}else{
  cvDF <- read_rds(paste0(outPrefix, "cvDF.rds"))
}

# gi <- read_rds(paste0(outPrefix, ".gi.rds"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))
# tidyCV <- read_feather(paste0(outPrefix, ".tidyCV.feather"))
# cvDF <- read_rds(paste0(outPrefix, "cvDF.rds"))
# designDF <- read_rds(paste0(outPrefix, "designDF.rds"))
#===============================================================================
# Performance Evaluation
#===============================================================================

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = cvDF$pred,
  labels = cvDF$label,
  modnames = str_c(cvDF$name, "_on_", cvDF$dep_var),
  dsids = cvDF$id,
  posclass = levels(cvDF$label[[1]])[2],
  x_bins = 100)

write_rds(curves, paste0(outPrefix, ".curves.rds"))
# curves <- read_rds(paste0(outPrefix, ".curves.rds"))

# get data.frame with auc values
aucDF <- as_tibble(auc(curves)) %>% 
  separate(modnames, into = c("name", "dep_var"), sep = "_on_") %>% 
  mutate(
    name = factor(name, unique(name)),
    dep_var = factor(str_replace(dep_var, "Loop_", ""), 
                    str_replace(unique(dep_var), "Loop_", ""))
         ) %>% 
  arrange(dep_var, name)

aucDFmed <- aucDF %>%
  group_by(dep_var, name, curvetypes) %>% 
  summarize(
    aucs_median = median(aucs, na.rm = TRUE),
    aucs_mean = mean(aucs, na.rm = TRUE),
    aucs_sd = sd(aucs, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(
    curvetypes = factor(curvetypes, c("ROC", "PRC"))
  )

write_feather(aucDFmed, paste0(outPrefix, ".aucDFmed.feather"))
# aucDFmed <- read_feather(paste0(outPrefix, ".aucDFmed.feather"))

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = name, y = aucs_mean, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) +
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, vjust = 1.5) +
  facet_grid(curvetypes ~ dep_var, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.barplot.pdf"), w = 14, h = 7)


# Only AUC PRC as barplot
p <- ggplot(filter(aucDFmed, curvetypes == "PRC"), 
            aes(x = name, y = aucs_mean, fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  facet_grid(. ~ dep_var) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance\n(auPRC)")

ggsave(p, file = paste0(outPrefix, ".AUC_PRC.barplot.pdf"), w = 14, h = 5)


#===============================================================================
# Analyse differnt subset of negativ loops depending on CTCF binding -----------
#===============================================================================
# There are different cases of non-looping pairs of CTCF sites: 
# - some are genomic sequences which are not even bound by CTCF, 
# - others are bound by CTCF but not engaged in any loop, 
# - and others are engaged in some loops, but do not form particular loop.

# Split CTCF pairs in the following classes
# - not CTCF bound (both?)
# - bound but not in any loop (both?)
# - bound and in some other loops (both?)

# get all anchor regions
reg <- regions(gi)

# get subset of interactions that are true loops
loops <- gi[gi$loop == "Loop"]

# get index of regions that are involve in true loops
loop_idx1 <- anchorIds(loops, "first")
loop_idx2 <- anchorIds(loops, "second")

# get subset of anchor regions that are involved in any true loop
loops_anchors <- c(reg[loop_idx1], reg[loop_idx2]) %>% unique()

# get boolean vector of interactions if they contain an anchor that is invovled in any peak
anchor_in_loop <- overlapsAny(gi, loops_anchors)



# load CTCF peak data as GRanges object
peak_df <- read_tsv((CTCF_peaks_bed), col_names = FALSE)
peak_gr <- GRanges(peak_df$X1, IRanges(peak_df$X2, peak_df$X3), strand = "*", 
                seqinfo = seqinfo(gi))

# get boolean vector of interactions if they overlap any peak
anchor_in_peak <- overlapsAny(gi, peak_gr)

# add columns to interaction df
df <- df %>% 
  mutate(
    anchor_in_peak = anchor_in_peak,
    anchor_in_loop = anchor_in_loop,
    
    # consturct negative set
    # - some are genomic sequences which are not even bound by CTCF, 
    loop_not_bound = loop == "Loop" | !anchor_in_peak,
    loop_not_bound = ifelse(loop_not_bound,
                            loop == "Loop",
                            NA) %>% 
      factor(c(FALSE, TRUE), LOOP_LEVEL), 
    
    # - others are bound by CTCF but not engaged in any loop, 
    loop_bound_and_not_anchor = loop == "Loop" | (anchor_in_peak & !anchor_in_loop),
    loop_bound_and_not_anchor = ifelse(loop_bound_and_not_anchor,
                            loop == "Loop",
                            NA) %>% 
      factor(c(FALSE, TRUE), LOOP_LEVEL), 
    
    # - and others are engaged in some loops, but do not form particular loop.
    loop_anchor_in_some_loops = loop == "Loop" | anchor_in_loop,
    loop_anchor_in_some_loops = ifelse(loop_anchor_in_some_loops,
                                       loop == "Loop",
                                       NA) %>% 
      factor(c(FALSE, TRUE), LOOP_LEVEL), 
  )

#===============================================================================
# Cross validation
#===============================================================================

# k fold cross validation with 1 repeats
set.seed(3579)

tidyCV <- read_feather(paste0(outPrefixSelectedModel, ".tidyCV.feather"))

dependent_variable <-  c("loop", "loop_not_bound", "loop_bound_and_not_anchor", 
                         "loop_anchor_in_some_loops")

designDF <- tibble(
  dep_var = dependent_variable,
  design = map(dep_var, function(dep_var) 
    c(
      list(
        "Dist" =  as.formula(str_c(dep_var, " ~ dist")),
        "Orientation" =  as.formula(str_c(dep_var, " ~ strandOrientation")),
        "Motif" =  as.formula(str_c(dep_var, " ~ score_min")),
        "Dist+Orientation+Motif" =  as.formula(str_c(dep_var, " ~ dist + strandOrientation + score_min"))
      ),
      map(meta$name, ~as.formula(str_c(dep_var, " ~ dist + strandOrientation + score_min + cor_", .x)) )
    )
  ),
  name = list(c("Dist", "Orientation", "Motif", "Dist+Orientation+Motif", meta$name)),
  color = list(c("greenyellow", "gold2", "khaki", "brown4", COL_SELECTED_TF_2))
) %>% 
  unnest(design, name, color) %>% 
  mutate(
    dep_var = factor(dep_var, unique(dep_var)),
    name = factor(name, unique(name))
  ) %>% 
  arrange(dep_var, name)

write_rds(designDF, paste0(outPrefix, "anchor_type_designDF.rds"))

# expand data.frame to have all combinations of model and split
cvDF <- tidyCV %>% 
  distinct(Fold) %>% 
  tidyr::expand(name = designDF$name, Fold) %>% 
  # add design formular for each TF
  left_join(designDF, by = "name") %>% 
  mutate(id = parse_integer(str_replace(Fold, "^Fold", "")))

if (!GI_LOCAL){
  # fit model on training part
  cvDF <- cvDF %>% 
    group_by(name, Fold) %>% 
    # fit model and save estimates in tidy format
    mutate(
      tidy_model = future_map2(Fold, design, .f = tidyer_fitter, 
                               tidyCV = tidyCV, data = df),
      pred = future_pmap(
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
  
  write_rds(cvDF, path = paste0(outPrefix, "anchor_type_cvDF.rds"))
  
}else{
  cvDF <- read_rds(paste0(outPrefix, "anchor_type_cvDF.rds"))
}

#===============================================================================
# Performance Evaluation
#===============================================================================

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = cvDF$pred,
  labels = cvDF$label,
  modnames = str_c(cvDF$name, "_on_", cvDF$dep_var),
  dsids = cvDF$id,
  posclass = levels(cvDF$label[[1]])[2],
  x_bins = 100)

write_rds(curves, paste0(outPrefix, ".anchor_type_curves.rds"))
# curves <- read_rds(paste0(outPrefix, ".anchor_type_curves.rds"))

# get data.frame with auc values
aucDF <- as_tibble(auc(curves)) %>% 
  separate(modnames, into = c("name", "dep_var"), sep = "_on_") %>% 
  mutate(
    name = factor(name, unique(name)),
    dep_var = factor(str_replace(dep_var, "loop_", ""), 
                     str_replace(unique(dep_var), "loop_", ""))
  ) %>% 
  arrange(dep_var, name)

aucDFmed <- aucDF %>%
  group_by(dep_var, name, curvetypes) %>% 
  summarize(
    aucs_median = median(aucs, na.rm = TRUE),
    aucs_mean = mean(aucs, na.rm = TRUE),
    aucs_sd = sd(aucs, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(
    curvetypes = factor(curvetypes, c("ROC", "PRC"))
  )

write_feather(aucDFmed, paste0(outPrefix, ".anchor_type_aucDFmed.feather"))
# aucDFmed <- read_feather(paste0(outPrefix, ".anchor_type_aucDFmed.feather"))

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDFmed, aes(x = name, y = aucs_mean, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) +
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  # geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), size = 3, vjust = 1.5) +
  facet_grid(curvetypes ~ dep_var, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".anchor_type_AUC_ROC_PRC.barplot.pdf"), w = 14, h = 7)


# Only AUC PRC as barplot
p <- ggplot(filter(aucDFmed, curvetypes == "PRC"), 
            aes(x = name, y = aucs_mean, fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = aucs_mean - aucs_sd, ymax = aucs_mean + aucs_sd),
                width = .25, position = position_dodge(width = 1)) + 
  geom_text(aes(label = round(aucs_mean, 2), y = aucs_mean - aucs_sd), 
            size = 5, hjust = 1.2, angle = 90) +
  facet_grid(. ~ dep_var) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = designDF$color) +
  labs(x = "Models", y = "Prediction performance\n(auPRC)")

ggsave(p, file = paste0(outPrefix, ".anchor_type_AUC_PRC.barplot.pdf"), w = 14, h = 5)

