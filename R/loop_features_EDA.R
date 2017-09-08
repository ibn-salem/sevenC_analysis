################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


# require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(RColorBrewer)   # for nice colors
require(feather)      # for efficient storing of data.frames
require(InteractionSet)


# 0) Set parameter --------------------------------------------------------

# use previously saved gi object?

MIN_MOTIF_SIG <- 6
WINDOW_SIZE <- 1000
BIN_SIZE <- 1


COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

outPrefix <- file.path("results", paste0("v04_selected_models.", 
                                         paste0("motifSig", MIN_MOTIF_SIG), 
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

COL_ANCHOR = brewer.pal(12, "Paired")[c(2,6)] #[c(2,1,5,6)]

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
  select(TF, name, filePath, everything())


#-------------------------------------------------------------------------------
# load main data set with loops and features
#-------------------------------------------------------------------------------
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

#===============================================================================
# Plot ChIP-seq coverage at anchors
#===============================================================================
# get list of coverage for all anchors

FACT = "ZNF143"
N_EXAMPLE = 3
covList <- mcols(regions(gi))[,paste0("cov_", FACT)]

# get randomly three loop and three non-loping pairs
bothSupport <- gi$Loop_Rao_GM12878  == "Loop" & gi$Loop_Tang2015_GM12878 == "Loop"
bothNoSupport <- gi$Loop_Rao_GM12878  != "Loop" & gi$Loop_Tang2015_GM12878 != "Loop"

rand_examples = c(
  sample(which(bothSupport), N_EXAMPLE),
  sample(which(bothNoSupport), N_EXAMPLE)
)
write_rds(rand_examples, paste0(outPrefix, ".rand_examples.rds"))

# Loop examples : 357242
# No Loop example : 214045

cov_up <- covList[anchors(gi[rand_examples], type = "first", id = TRUE)]
cov_down <- covList[anchors(gi[rand_examples], type = "second", id = TRUE)]

covDF <- tibble(
    interaction_id = rand_examples,
    loop = rep(c("Loop", "No loop"), each = 3),
    reg = c(1:N_EXAMPLE, 1:N_EXAMPLE),
    cov_up = as.list(cov_up),
    cov_down = as.list(cov_down),
  ) %>% 
  gather(key = "anchor", value = "cov", cov_up, cov_down) %>% 
  mutate(anchor = str_replace(anchor, "cov_", "")) %>% 
  unnest(cov) %>% 
  mutate(pos = rep(1:WINDOW_SIZE, n()/WINDOW_SIZE) - WINDOW_SIZE/2)

write_feather(covDF, paste0(outPrefix, ".6_random_pairs.feather"))

#-------------------------------------------------------------------
# plotting coverage 
#-------------------------------------------------------------------
p = ggplot(covDF, aes(x = pos, fill = anchor)) +
  geom_ribbon(aes(ymin = 0, ymax = cov)) +
  facet_grid(reg ~ loop, scales = "free_y") +
  theme_bw() + 
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text.y = element_blank()) + 
  labs(y = paste(FACT, "ChIP-seq"), x = "Positions around CTCF motif") +
  scale_fill_manual(values = alpha(COL_ANCHOR, 0.5)) 
# + 
  # coord_cartesian(ylim = c(0,8*mean(factDF$signal)))
# p
ggsave(p, file = paste0(outPrefix, ".6_random_pairs.pdf"), w = 7, h = 7)

#-------------------------------------------------------------------
# plot only two pairs sites
#-------------------------------------------------------------------
p = ggplot(filter(covDF, reg == 1), aes(x = pos, fill = anchor)) +
  geom_ribbon(aes(ymin = 0, ymax = cov)) +
  facet_grid(reg ~ loop, scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_blank()) + 
  labs(y = paste(FACT, "ChIP-seq"), x = "Positions around CTCF motif") +
  scale_fill_manual(values = alpha(COL_ANCHOR, 0.5)) 
ggsave(p, file = paste0(outPrefix, ".2_random_pairs.pdf"), w = 7, h = 3.5)



#-------------------------------------------------------------------
# factAncDF <- cast(sub3, TF+pos+reg+support~anchor, value="signal")
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
  theme_bw() + theme(strip.text.y = element_blank(), text = element_text(size = 20)) + 
  xlab("UP anchor coverage") + ylab("DOWN anchor coverage") + 
  # scale_color_gradient2(low = "grey", mid = scales::muted("blue"), high = "gray")
  # scale_color_gradientn(colors = RColorBrewer::brewer.pal(3, "BrBG"))
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
  theme_bw() + theme(strip.text.y = element_blank(), text = element_text(size = 20)) + 
  xlab("UP anchor coverage") + ylab("DOWN anchor coverage") + 
  # scale_color_gradient2(low = "grey", mid = scales::muted("blue"), high = "gray")
  # scale_color_gradientn(colors = RColorBrewer::brewer.pal(3, "BrBG"))
  scale_color_gradientn(colors = colorspace::rainbow_hcl(20))

ggsave(p, file = paste0(outPrefix, ".2_random_pairs.cor.pdf"), w = 7, h = 3.5)


#===============================================================================
# Compare Correlation with Boxplot
#===============================================================================

p <- ggplot(tidyDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = interaction(loop, name))) + 
  geom_boxplot(fill = "white", width = .2) +
  facet_grid(. ~ name) + 
  scale_fill_manual(values = COL_SELECTED_TF) +
  theme_bw() + 
  theme(
    text = element_text(size=20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "ChIP-seq\n Correlation", x="")
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



