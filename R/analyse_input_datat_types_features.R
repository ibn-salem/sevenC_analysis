#*******************************************************************************
# Analysis different input data types for loop prediction with chromloop. ------
#*******************************************************************************
library(tidyverse)    # for tidy data
library(RColorBrewer)   # for nice colors
library(feather)      # for efficient storing of data.frames

# Set parameter ----------------------------------------------------------------

# use previously saved gi object?
GI_LOCAL <- FALSE
N_CORES = min(10, parallel::detectCores() - 1)

# MIN_MOTIF_SIG <- 5
MOTIF_PVAL <- 2.5 * 1e-06
WINDOW_SIZE <- 1000
BIN_SIZE <- 1
K = 10  # K-fold corss validation
N_TOP_MODELS = 10

outPrefix <- file.path("results", paste0("v05_input_types.", 
                                         paste0("motifPval", MOTIF_PVAL), 
                                         "_w", WINDOW_SIZE, 
                                         "_b", BIN_SIZE))

COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")


# Parse and prcessed data and meta data  ---------------------------------------
meta <- read_tsv(paste0(outPrefix, ".meta_filtered.tsv"))

df <- read_feather(paste0(outPrefix, ".df.feather"))


# make a tidy DF
tidyDF <- df %>% 
  gather(starts_with("cor_"), key = name, value = cor) 
  # mutate(type = str_replace(type, "^cor_", ""))
  mutate(name = str_sub(name, 5))

# Compare Correlation with Boxplot ---------------------------------------------

tidySubDF <- tidyDF %>% 
  sample_n(min(nrow(tidyDF), 10^7))

p <- ggplot(tidySubDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = name), lwd = 1.5) + 
  geom_boxplot(fill = "white", lwd = 1.5, width = .2) +
  facet_grid(. ~ name) + 
  scale_fill_manual(
    values = colorRampPalette(brewer.pal(12, "Set3"))(nrow(meta)), 
    guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size = 20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Correlation of ChIP-seq signal", x = "")
# geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,3)), x=1.5, y=1.1), size=5)
# p
ggsave(p, file = paste0(outPrefix, ".cor.by_name_and_loop.boxplot.pdf"), w = 14, h = 7)


#======================== Compare Correlation with Boxplot ==================

p <- ggplot(tidyDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = interaction(loop, name))) + 
  geom_boxplot(fill = "white", width = .2) +
  facet_grid(. ~ name) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(nrow(meta))) +
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
    text = element_text(size = 20),
    # axis.text.x=element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "ChIP-seq\n Correlation", x = "")
ggsave(p, file = paste0(outPrefix, ".EDA.cor.by_TF_and_loop.boxplot.pdf"), w = 7, h = 3.5)

