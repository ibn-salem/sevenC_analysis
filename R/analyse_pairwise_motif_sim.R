#*******************************************************************************
# Anaylse pairwise motif similairity by alignemnt for coandiate pairs as feature
#*******************************************************************************

library(sevenC)  # to import() BED files
require(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()
require(RColorBrewer)   # for nice colors
library(ggsignif)

#*******************************************************************************
# Parameters ----
#*******************************************************************************
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

MOTIF_PVAL <- 2.5 * 1e-06

SHAPE_WINDOW = 500
SEQ_WINDOW = 19

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", paste0("CTCF_JASPAR.v01.motifSeqSim.sw", SEQ_WINDOW))


#*******************************************************************************
# Analze moitf seq similarity
#*******************************************************************************
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi.rds"))

BSgenome = BSgenome.Hsapiens.UCSC.hg19



getPairwiseDist <- function(gi, BSgenome, window_size){

  # extract sequences for all anchors
  windowGR <- resize(regions(gi), width = window_size, fix = "center")
  ancSeq <- getSeq(BSgenome, windowGR)
  
  editDist <- map2_dbl(
    .x = anchors(gi, "first", id = TRUE),
    .y = anchors(gi, "second", id = TRUE),
    .f = ~ stringDist(c(ancSeq[.x], ancSeq[.y]))[1]
  )
  return(editDist)
}

getPairwisePID <- function(gi, BSgenome, window_size){
  
  # extract sequences for all anchors
  windowGR <- resize(regions(gi), width = window_size, fix = "center")
  ancSeq <- getSeq(BSgenome, windowGR)
  
  identity <- map2_dbl(
    .x = anchors(gi, "first", id = TRUE),
    .y = anchors(gi, "second", id = TRUE),
    .f = ~ pid(pairwiseAlignment(ancSeq[.x], ancSeq[.y]))
  )
  return(identity)
}

getAlignScore <- function(gi, BSgenome, window_size){
  
  # extract sequences for all anchors
  windowGR <- resize(regions(gi), width = window_size, fix = "center")
  ancSeq <- getSeq(BSgenome, windowGR)
  
  score <- map2_dbl(
    .x = anchors(gi, "first", id = TRUE),
    .y = anchors(gi, "second", id = TRUE),
    .f = ~ pairwiseAlignment(ancSeq[.x], ancSeq[.y], scoreOnly = TRUE)
  )
  return(identity)
}


n = 20000
subGI <- gi[sample.int(length(gi), size = n)]

# add edit distance
subGI$seqDist <- getPairwiseDist(subGI, BSgenome, window_size = SEQ_WINDOW)
subGI$pid <- getPairwisePID(subGI, BSgenome, window_size = SEQ_WINDOW)
subGI$alignScore <- getAlignScore(subGI, BSgenome, window_size = SEQ_WINDOW)

p <- ggplot(
  as.data.frame(mcols(subGI)), aes(y = seqDist, x = loop, color = loop)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Loop", "No loop")), color = "black") +
  ylim(0, 15) +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Sequence\ndistance", x = "")
ggsave(p, file = paste0(outPrefix, ".", n, ".editDistance_by_loop.boxplot.pdf"), w = 3, h = 3)

p <- ggplot(
  as.data.frame(mcols(subGI)), aes(y = alignScore, x = loop, color = loop)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Loop", "No loop")), color = "black") +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "Alignment score", x = "")
ggsave(p, file = paste0(outPrefix, ".", n, ".alignScore_by_loop.boxplot.pdf"), w = 3, h = 3)


countDF <- as.data.frame(mcols(subGI)) %>% 
  mutate(seqDist = as.integer(seqDist)) %>% 
  group_by(loop, seqDist) %>% 
  summarize(
    n = n()
  ) %>% 
  mutate(percent = n / sum(n) * 100)

p <- ggplot(countDF, aes(y = percent, x = seqDist, fill = loop)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size = 20),
    legend.position = "bottom") + 
  labs(y = "%", x = "Sequence distance")

ggsave(p, file = paste0(outPrefix, ".", n, ".editDistance_by_loop.barplot.pdf"), w = 3, h = 3)


pdf(paste0(outPrefix, ".", n, ".editDistance_by_loop.pdf"), w = 3, h = 6)
  boxplot(seqDist ~ loop, data = mcols(subGI))
dev.off()


