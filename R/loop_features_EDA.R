################################################################################
# Analysis of predictd chromatin looping interactions using the chromloop tool
################################################################################


# require(chromloop)    # devtools::install_github("ibn-salem/chromloop")
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(RColorBrewer)   # for nice colors
require(feather)      # for efficient storing of data.frames


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

#===============================================================================
# Plot ChIP-seq coverage at anchors
#===============================================================================


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



