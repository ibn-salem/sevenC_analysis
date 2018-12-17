#*******************************************************************************
# Compare the performance of sevenC to other tools
#*******************************************************************************

library(sevenC)    # devtools::install_github("ibn-salem/sevenC")
library(tidyverse)    # for tidy data
library(precrec)      # for ROC and PRC curves
library(RColorBrewer)   # for nice colors
library(feather)      # for efficient storing of data.frames
library(ROCR)         # for binary clasification metrices

source("R/sevenC.functions.R")

# Set parameter ----------------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(16, parallel::detectCores() - 1)

# MIN_MOTIF_SIG <- 5
MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

Oti2016_file <- "data/ENCODE/Peaks/ENCFF710VEH.bed.CTCF_motif.bed.Oti_peaks.bed"

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

out_prefix_selected_models <- file.path("results", paste0("v05_selected_models.", 
                                                          paste0("motifPval", MOTIF_PVAL), 
                                                          "_w", WINDOW_SIZE, 
                                                          "_b", BIN_SIZE))

screenPrefix <- file.path("results", paste0("v05_screen_TF_lfc.", 
                                            paste0("motifPval", MOTIF_PVAL), 
                                            "_w", WINDOW_SIZE, 
                                            "_b", BIN_SIZE))

outPrefix <- file.path("results", paste0("v05_compare_tools.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

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

# Select motifs and parse input data -----------------------------------

# read CTCF moitf pairs as candidates
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))

# rad Oti2016 loops -------------------------------------

read_Oti2016 <- function(in_file){
  
  df <- read_tsv(Oti2016_file, skip = 1, col_names = FALSE) %>% 
    mutate(
      up_anc = X4 %>% str_split_fixed("-", 2) %>% magrittr::extract(, 1),
      down_anc = X4 %>% str_split_fixed("-", 2) %>% magrittr::extract(, 2),
      up_chr = up_anc %>% str_extract("chr[^:]*"),
      up_start = up_anc %>% str_match(":(\\d*)\\.\\.") %>% magrittr::extract(, 2) %>% parse_integer(),
      up_end = up_anc %>% str_match("\\.\\.(\\d*)") %>% magrittr::extract(, 2) %>% parse_integer(),
      down_chr = down_anc %>% str_extract("chr[^:]*"),
      down_start = down_anc %>% str_match(":(\\d*)\\.\\.") %>% magrittr::extract(, 2) %>% parse_integer(),
      down_end = down_anc %>% str_match("\\.\\.(\\d*)") %>% magrittr::extract(, 2) %>% parse_integer()
    )
  
  gr1 <- GRanges(df$up_chr, IRanges(df$up_start, df$up_end))
  gr2 <- GRanges(df$down_chr, IRanges(df$down_start, df$down_end))
  gi <- GInteractions(gr1, gr2, score = df$X5)
  
}

Oti_loops <- read_Oti2016(Oti2016_file)

mcols(gi)$oti <- overlapsAny(gi, Oti_loops)

# add score for each overlapping loop
ol <- findOverlaps(gi, Oti_loops)
mcols(gi)$Oti_score <- NA
mcols(gi)$Oti_score[queryHits(ol)] <- Oti_loops$score[subjectHits(ol)]


# add overlap with Oti loops
df <- mcols(gi) %>% 
  as.data.frame() %>%  
  as.tibble() 

write_rds(df, str_c(outPrefix, ".df.rds"))

# count predicted and true loops
df %>% dplyr::count(loop, oti) %>% 
  mutate(percent = n / sum(n) * 100)

# sens <- 9239 / (20786 + 9239) * 100
# PPV <- 9239 / (9239 + 5514) * 100

# Compare to 7C ----------------------------------------------------------------
cvDF <- read_rds(paste0(out_prefix_selected_models, "cvDF.rds"))

cvDF <- cvDF %>%
  filter(!str_detect(name, ".*_only$") )

# combine folds from cross-validation
predDF <- cvDF %>% 
  group_by(name) %>% 
  unnest(pred, label) %>% 
  group_by(name, color) %>% 
  summarize(pred = list(pred), label = list(label))

# add Oti data
predDF <- predDF %>% 
  ungroup() %>% 
  mutate(name = as.character(name)) %>% 
  bind_rows(
    tibble(name = "Oti2016", color = RColorBrewer::brewer.pal(8, "Set1")[8],
           pred = list(df$Oti_score), label = list(df$loop))
  )

# get AUC of ROC and PRC curves for all 
curves <- evalmod(
  scores = predDF$pred,
  labels = predDF$label,
  modnames = predDF$name,
  posclass = levels(predDF$label[[1]])[2],
  x_bins = 100
)

#autoplot(curves)
write_rds(curves, paste0(outPrefix, ".curves.rds"))
# curves <- read_rds(paste0(outPrefix, ".curves.rds"))

# get data.frame with auc values
aucDF <-  as_tibble(auc(curves)) %>% 
  mutate(
    modnames = factor(modnames, unique(modnames)),
    curvetypes = factor(curvetypes, c("ROC", "PRC"))
    ) %>% 
  arrange(modnames)

#-------------------------------------------------------------------------------
# barplot of AUCs of ROC and PRC
#-------------------------------------------------------------------------------
p <- ggplot(aucDF, aes(x = modnames, y = aucs, fill = modnames, label = round(aucs, 2))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(size = 5, hjust = "top", angle = 90) +
  facet_grid(curvetypes ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = predDF$color) +
  labs(x = "Models", y = "Prediction performance (AUC)")
ggsave(p, file = paste0(outPrefix, ".AUC_ROC_PRC.barplot.pdf"), w = 5, h = 7)

#-------------------------------------------------------------------------------
# ROC and PRC curves
#-------------------------------------------------------------------------------
curveDF <- as_tibble(as.data.frame(curves)) %>% 
  mutate(modname = factor(modname, unique(modname))) %>% 
  # downsmaple points from ROC curve for to fix large plot file size
  sample_n(20000)

# curveDF <- 

roc_labels <- aucDF %>% 
  filter(curvetypes == "ROC") %>% 
  mutate(label = str_c(modnames, " (auROC=", signif(aucs, 2), ")")) %>% 
  pull(label)
         
  
p <- ggplot(filter(curveDF, type == "ROC"), aes(x = x, y = y, color = modname)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dotted") +
  theme_bw() + theme(aspect.ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual("", values = predDF$color, labels = roc_labels) +
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 8),
        legend.position = "right")

ggsave(p, file = paste0(outPrefix, ".ROC_with_legned.pdf"), w = 5, h = 5)
ggsave(p + theme(legend.position = "none"), file = paste0(outPrefix, ".ROC.pdf"), w = 5, h = 5)


# get PRC plots

prc_base = precrec:::.get_pn_info(curves)$prc_base

prc_labels <- aucDF %>% 
  filter(curvetypes == "PRC") %>% 
  mutate(label = str_c(modnames, " (auPRC=", signif(aucs, 2), ")")) %>% 
  pull(label)

p <- ggplot(filter(curveDF, type == "PRC"), aes(x = x, y = y, color = modname)) +
  geom_line() +
  geom_abline(intercept = prc_base, slope = 0, lty = "dotted") +
  theme_bw() + theme(aspect.ratio = 1) +
  labs(x = "Recall", y = "Precision") +
  scale_color_manual("", values = predDF$color, labels = roc_labels) + 
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 8),
        legend.position = "right")

ggsave(p, file = paste0(outPrefix, ".PRC_with_legned.pdf"), w = 5, h = 5)
ggsave(p + theme(legend.position = "none"), file = paste0(outPrefix, ".PRC.pdf"), w = 5, h = 5)

# Binary predictions ===========================================================

# compute binary prediction by appling 7C default cutoff and 0 for Oti2016

binaryDF <- predDF %>% 
  filter(name %in% c(SELECTED_TF, "Oti2016")) %>% 
  mutate(cutoff = c(rep(sevenC::cutoffBest10, length(SELECTED_TF)), 0)) %>% 
  mutate(
    name = factor(ifelse(name != "Oti2016", str_c("7C_", name), name), 
                  c(str_c("7C_", SELECTED_TF), "Oti2016")),
    pred_label = map2(pred, cutoff, 
                      ~factor(!is.na(.x) & .x >= .y,
                                   c(FALSE, TRUE), 
                                   c("No loop", "Loop")))) %>% 
  select(name, color, pred_label, label) %>% 
  unnest(pred_label, label)

# compute binary classification evaluation performance
manual_performance <- binaryDF %>% 
  mutate(pred = pred_label == "Loop",
         true = label == "Loop") %>% 
  group_by(name) %>% 
  summarize(
    TP = sum(pred & true),
    FP = sum(pred & !true),
    TN = sum(!pred & !true),
    FN = sum(!pred & true),
    N = TP + FP + TN + FN,
    sensitivity = TP / (TP + FN) * 100,
    specificity = TN / (TN + FP) * 100,
    precision = TP / (TP + FP) * 100,
    accuracy = (TP + TN) / N * 100
  )
  
write_tsv(manual_performance, str_c(outPrefix, ".binary.manual_performance.tsv"))

# plot binary prediction performance
p <- manual_performance %>% 
  gather("metric", "performance", sensitivity, specificity, precision, accuracy) %>% 
  ggplot(aes(x = name, y = performance, fill = name, 
                               label = round(performance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = "top") +
  facet_grid(metric ~ ., scale = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c(unique(binaryDF$color), brewer.pal(8, "Set1")[8])) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "")

ggsave(p, file = str_c(outPrefix, ".binary_performance.barplot.pdf"), w = 5, h = 5)
