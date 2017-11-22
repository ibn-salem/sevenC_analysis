# Analysis of predictd chromatin looping interactions using the chromloop tool


# require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(InteractionSet)
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(RColorBrewer)   # for nice colors
require(feather)      # for efficient storing of data.frames


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?

# MIN_MOTIF_SIG <- 5
MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1

dataCandidatesPreifx <- file.path("results", 
                                 paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

# outPrefix <- file.path("results", paste0("v05_selected_models.", 
#                                          paste0("motifSig", MIN_MOTIF_SIG), 
#                                          "_w", WINDOW_SIZE, 
#                                          "_b", BIN_SIZE))
outPrefix <- file.path("results", paste0("v05_selected_models.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

dir.create(dirname(outPrefix), showWarnings = FALSE)

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

# COL_SELECTED_TF = brewer.pal(length(SELECTED_TF), "Set1")
COL_SELECTED_TF = brewer.pal(12, "Paired")
COL_SELECTED_TF_1 = brewer.pal(12, "Paired")[c(1, 3, 5, 7, 9, 11)]
COL_SELECTED_TF_2 = brewer.pal(12, "Paired")[c(2, 4, 6, 8, 10, 12)]

COL_ANCHOR = brewer.pal(12, "Paired")[c(6,2)] #[c(2,6)] #[c(2,1,5,6)]

#-------------------------------Parse and filter input ChiP-seq data  

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
  select(TF, name, filePath, everything())


# ----------------load main data set with loops and features
# write_feather(df, paste0(outPrefix, ".df.feather"))
# df <- read_feather(paste0(outPrefix, ".df.feather"))
# write_feather(edaDF, paste0(outPrefix, ".edaDF.feather"))

df <- read_feather(paste0(outPrefix, ".df.feather"))

# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = name, value = cor) %>% 
  mutate(name = str_sub(name, 5)) %>% 
  mutate(name = factor(name, SELECTED_TF))

# tidySeltDF <- tidyDF %>% 
#   filter(name %in% SELECTED_TF)

# load GInteractions object
gi <- read_rds(paste0(outPrefix, ".gi.rds"))



#*******************************************************************************
# Compaire motifs from RSAT and JASPAR ----
#*******************************************************************************
jasparGR <- read_rds(paste0(dataCandidatesPreifx, "motifGR.rds")) 
rsatGR <- chromloop::motif.hg19.CTCF

jaspar_unique <- sum(countOverlaps(jasparGR, rsatGR) == 0)
rsat_unique <- sum(countOverlaps(rsatGR, jasparGR) == 0)

jaspar_common <- sum(countOverlaps(jasparGR, rsatGR) > 0)
rsat_common <- sum(countOverlaps(rsatGR, jasparGR) > 0)

motif_counts <- tibble(
  type = c("JASPAR only", "Common", "RSAT only"),
  count = c(jaspar_unique, jaspar_common, rsat_unique)
)
p <- ggplot(motif_counts, aes(x = type, y = count)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = count), vjust = "bottom")
ggsave(paste0(outPrefix, ".motif_overlap_JASPAR_RSAT.barplot.pdf"), w = 3, h = 3)
# -------------------

#======================== Percent positives  ==================

# simple stats:
countsDF <- tibble(
  nMotif = length(regions(gi)),
  nPairs = nrow(df),
  HiC = sum(df$Loop_Rao_GM12878 == "Loop"),
  HiC_percent = HiC / nPairs * 100, 
  HIC_ChIAPET = sum(df$Loop_Tang2015_GM12878 == "Loop" | df$Loop_Rao_GM12878 == "Loop"),
  HIC_ChIAPET_percent = HIC_ChIAPET / nPairs * 100,
  # HIC_ChIAPET_CaptureC = sum(df$HIC_ChIAPET_CaptureC == "Loop"),
  # HIC_ChIAPET_CaptureC_percent = HIC_ChIAPET_CaptureC / nPairs * 100
  loop = sum(df$loop == "Loop"),
  loop_percent = loop / nPairs * 100,
)

write_tsv(countsDF, path = paste0(outPrefix, ".positives.countsDF.tsv"))

# vennList <- map(
#   select(df, starts_with("Loop_")),
#   function(x) which(x == "Loop")
# )


# #--------- VennDiagram of loops -----------------------------------
# require(VennDiagram)# for VennDiagrams
# venn.diagram(
#   x = vennList, 
#   filename = paste0(outPrefix, ".loop_balance.venn.png"),
#   imagetype = "png",
#   euler.d = TRUE, scaled = TRUE
# )
# 
# require(venneuler)  # for Venn-Euler diagram (area propotional)
# 
# vennMat <- df %>% 
#   select(starts_with("Loop_")) %>% 
#   mutate_all( function(x) x == "Loop")
# 
# pdf(paste0(outPrefix, ".loop_balance.venneuler.pdf"))
# plot(venneuler(vennMat))
# dev.off()


# Venn-Diagram of ChIA-PET and Hi-C loops --------------------------------------

# vennList <- map(
#   select(df, Loop_Rao_GM12878, Loop_Tang2015_GM12878_CTCF, Loop_Tang2015_GM12878_PolII),
#   function(x) which(x == "Loop")
# )
# names(vennList) <- c("Hi-C", "CTCF ChIA-PET", "PolII ChIA-PET")
# 
# venn.diagram(
#   x = vennList, 
#   filename = paste0(outPrefix, ".loop_balance.Hi-C_ChIA-PET.venn.png"),
#   imagetype = "png",
#   euler.d = TRUE, scaled = TRUE,
#   lwd = 6,
#   lty = 'blank',
#   fill = c("cornflowerblue", "green", "yellow"),
#   cex = 2,
#   fontface = "bold",
#   fontfamily = "sans",
#   cat.cex = 2.5,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   rotation = 1,
#   cat.just = list(c(0.6,1) , c(.7,1.2) , c(1,0))
# )

# vennMat <- df %>% 
#   select(Loop_Rao_GM12878, Loop_Tang2015_GM12878_CTCF, Loop_Tang2015_GM12878_PolII) %>% 
#   mutate_all( function(x) x == "Loop")
# names(vennMat) <- c("Hi-C", "CTCF ChIA-PET", "PolII ChIA-PET")
# 
# pdf(paste0(outPrefix, ".loop_balance.Hi-C_ChIA-PET.venneuler.pdf"))
# plot(venneuler(vennMat))
# dev.off()


# ---------------- Pie chart of percent positives ---------------------


nDF <- df %>% 
  group_by(loop) %>% 
  summarize(
    n = n()
    ) %>% 
  mutate(percent = n/sum(n) * 100)

p <- ggplot(nDF, aes( x = "", y = n, fill = loop)) +
  geom_bar(stat = "identity", width = 0.5) + 
  coord_polar("y") +
  scale_fill_manual(values = COL_LOOP) + 
  geom_text(aes(y = cumsum(rev(n) / 2),
                label = paste0(rev(loop), "\n", rev(n), "\n", 
                              rev(round(percent, 2)), "%")),
            size = 5) + 
  geom_text(aes(x = 0.5, y = 0, label = paste("n =", nrow(df))), size = 5) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )
p
ggsave(p, file = paste0(outPrefix, ".loop_balance.pie.pdf"), w = 6, h = 6)


#======================== distance of looping and non looping pairs  ======


# plotDF <- melt(loopDF, id=c("IDg1","IDg1", "dist"), measure.vars=LOOP_COL_ALL)

p <- ggplot(df, aes(dist/10^3, fill = loop, color = loop)) +
  geom_density(size = 2, alpha = 0.5) +
  # facet_grid(variable ~ .) + 
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = "bottom") + 
  labs(y = "Density", x = "Loop Distance [kb]") 
ggsave(p, file = paste0(outPrefix, ".distance_by_loop.density.pdf"), w = 7, h = 3.5)


#======================== Motif orientation ======================

oreientationDF <- df %>%
  group_by(loop, strandOrientation) %>%
  summarise(
    n = n()
  ) %>%
  mutate(percent = n / sum(n) * 100)

p <- ggplot(oreientationDF, 
            aes(x = strandOrientation, y = n, fill = strandOrientation)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = 2) +
  facet_grid(loop ~ . , scales = "free_y") + 
  theme_bw()

ggsave(paste0(outPrefix, "orientation.n_by_orientation_and_loop.barplot.pdf"),
       w = 6, h = 6)


p <- ggplot(oreientationDF, 
            aes(x = loop, y = percent, fill = loop, color = loop)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(percent, 2)), size = 5, vjust = "inward", color = "black") +
  facet_grid(. ~ strandOrientation , scales = "free_y") + 
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = "bottom",
                     axis.text.x = element_blank()) +
  labs(x = "", y = "CTCF motif pairs [%]")

ggsave(paste0(outPrefix, "orientation_by_loop.percent.barplot.pdf"),
       w = 6, h = 3.5)


#======================== Motif score by loop ==================================
p <- ggplot(df, aes(x = score_min, fill = loop, color = loop)) +
  geom_density(size = 2, alpha = 0.5) +
  # facet_grid(variable ~ .) + 
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = "bottom") + 
  labs(y = "Density", x = "Motif hit significance [-log10(p)]")
p
ggsave(p, file = paste0(outPrefix, ".motif_score_by_loop.density.pdf"), w = 7, h = 3.5)

ancLoop <- unlist(anchors(gi[df$loop == "Loop"], id = TRUE))
# motif_sig <- mcols(regions(gi))[, "sig"]
motif_sig <- mcols(regions(gi))[, "score"]
motifDF <- tibble(
  motif_sig = motif_sig,
  loop = ifelse(1:length(motif_sig) %in% ancLoop, "Loop", "No loop")
)

p <- ggplot(motifDF, aes(x = motif_sig, fill = loop, color = loop)) +
  geom_density(size = 2, alpha = 0.5) +
  # facet_grid(variable ~ .) + 
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size = 20), legend.position = "bottom") + 
  labs(y = "Density", x = "Motif hit significance [-log10(p)]")
ggsave(p, file = paste0(outPrefix, ".anchor_motif_score_by_loop.density.pdf"), w = 7, h = 3.5)

# ====================  Plot ChIP-seq coverage at anchors ======================
# get list of coverage for all anchors

FACT = "ZNF143"
N_EXAMPLE = 3
# covList <- mcols(regions(gi))[, paste0("cov_", FACT)]
covList <- mcols(regions(gi))[, FACT]

# get randomly three loop and three non-loping pairs
bothSupport <- gi$Loop_Rao_GM12878  == "Loop" & gi$Loop_Tang2015_GM12878 == "Loop"
bothNoSupport <- gi$Loop_Rao_GM12878  != "Loop" & gi$Loop_Tang2015_GM12878 != "Loop"

rand_examples = c(
  sample(which(bothSupport), N_EXAMPLE),
  sample(which(bothNoSupport), N_EXAMPLE)
)

# rand_examples = c(
#   379220, 357242, 16818,
#   163948, 69243, 214045, 
# )
write_rds(rand_examples, paste0(outPrefix, ".rand_examples.rds"))

# Loop examples : 357242, 16818, 452694, 379220
# No Loop example : 214045, 69243, 163948

cov_up <- covList[anchors(gi[rand_examples], type = "first", id = TRUE)]
cov_down <- covList[anchors(gi[rand_examples], type = "second", id = TRUE)]

covDF <- tibble(
    interaction_id = rand_examples,
    loop = rep(c("Loop", "No loop"), each = 3),
    reg = c(3, 2, 1, 3, 2, 1),
    cov_up = as.list(cov_up),
    cov_down = as.list(cov_down),
  ) %>%
  gather(key = "anchor", value = "cov", cov_up, cov_down) %>%
  mutate(anchor = factor(str_replace(anchor, "cov_", ""), c("up", "down"))) %>%
  unnest(cov) %>%
  mutate(pos = rep(1:WINDOW_SIZE, n()/WINDOW_SIZE) - WINDOW_SIZE/2)

write_feather(covDF, paste0(outPrefix, ".6_random_pairs.feather"))

#---------------plotting coverage ----------------------------------------------

p = ggplot(covDF, aes(x = pos, fill = anchor)) +
  geom_ribbon(aes(ymin = 0, ymax = cov)) +
  facet_grid(reg ~ loop, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_blank()) +
  labs(y = paste(FACT, "ChIP-seq"), x = "Positions around CTCF motif") +
  scale_fill_manual(values = alpha(COL_ANCHOR, 0.5))
ggsave(p, file = paste0(outPrefix, ".6_random_pairs.pdf"), w = 7, h = 7)

#---------------------plot only two pairs sites ------------------------
p = ggplot(filter(covDF, reg == 1), aes(x = pos, fill = anchor)) +
  geom_ribbon(aes(ymin = 0, ymax = cov)) +
  facet_grid(reg ~ loop, scales = "free_y") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text.y = element_blank()) + 
  labs(y = paste(FACT, "ChIP-seq"), x = "Positions around CTCF motif") +
  scale_fill_manual(values = alpha(COL_ANCHOR, 0.5)) 
ggsave(p, file = paste0(outPrefix, ".2_random_pairs.pdf"), w = 7, h = 3.5)


factAncDF <- covDF %>% 
  spread(key = anchor, value = cov)


corDF <- factAncDF %>% 
  group_by(reg, loop) %>% 
  summarize(
    R = cor(up, down, use = "complete"),
    p_val = cor.test(up, down)$p.value
  )

p = ggplot(factAncDF, aes(x = up, y = down, color = pos)) +
  geom_point() +
  facet_grid(reg ~ loop, scales = "free") +
  geom_text(aes(x = max(factAncDF$up, na.rm = TRUE), y = 1, 
                label = paste0("R=", round(R, 2))),
            data = corDF, color = "black",
            vjust = "inward", hjust = "inward", size = 7.5) +
  theme_bw() + 
  theme(strip.text.y = element_blank(), 
        text = element_text(size = 20),
        axis.title.x = element_text(face = "bold", color = COL_ANCHOR[1]),
        axis.text.x = element_text(face = "bold", color = COL_ANCHOR[1]),
        axis.title.y = element_text(face = "bold", color = COL_ANCHOR[2]),
        axis.text.y = element_text(face = "bold", color = COL_ANCHOR[2])
        ) + 
  xlab("UP anchor coverage") + ylab("DOWN anchor coverage") + 
  scale_color_gradientn(colors = colorspace::rainbow_hcl(20))
ggsave(p, file = paste0(outPrefix, ".6_random_pairs.cor.pdf"), w = 7, h = 7)

subAncDF <- filter(factAncDF, reg == 1)

p = ggplot(filter(subAncDF, reg == 1), aes(x = up, y = down, color = pos)) +
  geom_point() +
  facet_grid(reg ~ loop, scales = "free") +
  geom_text(aes(x = max(subAncDF$up, na.rm = TRUE), y = 1, 
                label = paste0("R=", round(R, 2))),
            data =filter(corDF, reg == 1), color = "black",
            vjust = "inward", hjust = "inward", size = 7.5) +
  theme_bw() + theme(strip.text.y = element_blank(), 
                     axis.title.x = element_text(face = "bold", color = COL_ANCHOR[1]),
                     axis.text.x = element_text(face = "bold", color = COL_ANCHOR[1]),
                     axis.title.y = element_text(face = "bold", color = COL_ANCHOR[2]),
                     axis.text.y = element_text(face = "bold", color = COL_ANCHOR[2]),
                     text = element_text(size = 20)) + 
  xlab("UP anchor coverage") + ylab("DOWN anchor coverage") + 
  # scale_color_gradient2(low = "grey", mid = scales::muted("blue"), high = "gray")
  # scale_color_gradientn(colors = RColorBrewer::brewer.pal(3, "BrBG"))
  scale_color_gradientn(colors = colorspace::rainbow_hcl(20))
ggsave(p, file = paste0(outPrefix, ".2_random_pairs.cor.pdf"), w = 7, h = 3.5)


#======================== Compare Correlation with Boxplot ==================

p <- ggplot(tidyDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = interaction(loop, name))) + 
  geom_boxplot(fill = "white", width = .2) +
  facet_grid(. ~ name) + 
  scale_fill_manual(values = COL_SELECTED_TF) +
  theme_bw() + 
  theme(
    text = element_text(size = 20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "ChIP-seq\n Correlation", x = "")
# geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
# p
ggsave(p, file = paste0(outPrefix, ".EDA.cor.by_TF_and_loop.violin.pdf"), w = 7, h = 3.5)


p <- ggplot(tidyDF, aes(x = name, y = cor, col = loop)) +
  geom_boxplot() +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size=20),
    # axis.text.x=element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "ChIP-seq\n Correlation", x="")
ggsave(p, file = paste0(outPrefix, ".EDA.cor.by_TF_and_loop.boxplot.pdf"), w = 7, h = 3.5)



#-------------------------------------------------------------------------------
# Compare dist vs. log10(dist) in prediction 
#-------------------------------------------------------------------------------
stop("Stop here. Needs to be reimplemented.")
formula_list <- list(
  dist = as.formula("loop ~ dist + strandOrientation + score_min + cor_CTCF_lfc"),
  log10dist = as.formula("loop ~ dist_log10 + strandOrientation + score_min + cor_CTCF_lfc")
)

mod_list <- formula_list %>% 
  map(glm, family = binomial(), data = df)

param_list <- mod_list %>% 
  map(tidy) %>% 
  map("estimate")

pred_list <- map2(formula_list, param_list, 
                  ~ chromloop::pred_logit(data = df, formula = .x, betas = .y))


mdat <- mmdata(pred_list, labels = list(df$loop, df$loop), modnames = c("dist", "log10_dist"))
curves <- evalmod(mdat)

p <- autoplot(curves)
ggsave(p, paste0(outPrefix, ".compare_dist_vs_log10dist.corCTCF_lfc.pdf"))

