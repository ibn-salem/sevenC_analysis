#*******************************************************************************
# Anaylse DNA shape at loop anchors
#*******************************************************************************


library(chromloop)  # to import() BED files
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
library(tidyverse)    # for tidy data
library(stringr)      # for string functions
library(readr)        # for write_rds()
require(RColorBrewer)   # for nice colors
require(DNAshapeR)
require(BSgenome.Hsapiens.UCSC.hg19)

#*******************************************************************************
# Parameters ----
#*******************************************************************************
COL_LOOP = brewer.pal(8, "Dark2")[c(8,5)] #[c(2,1,5,6)]
names(COL_LOOP) <- c("No loop", "Loop")

MOTIF_PVAL <- 2.5 * 1e-06

SHAPE_WINDOW = 500

JASPAR_HG19_CTCF <- "data/JASPAR2018/MA0139.1.tsv"

# define data candidate path
dataCandidatesPreifx <- file.path("results", 
                                  paste0("CTCF_JASPAR.v01.pval_", MOTIF_PVAL))

outPrefix <- file.path("results", paste0("CTCF_JASPAR.v01.DNA_shape.w", SHAPE_WINDOW))


#*******************************************************************************
# Analze DNA sahpe
#*******************************************************************************
gi <- read_rds(paste0(dataCandidatesPreifx, ".gi_with_shape.rds"))

BSgenome = BSgenome.Hsapiens.UCSC.hg19



# define temp file for sequences
# tempFile <- tempfile("tmp.fa")
tempFile <- paste0(outPrefix, ".seq.fa")
                      
# get sequence of anchors
getFasta(regions(gi), BSgenome, width = SHAPE_WINDOW, filename = tempFile)

# predict DNA shape
shapeTypes <- c("MGW", "HelT", "ProT", "Roll", "EP", "Opening", "Rise", "Shift", "Stagger", "Slide")
# shapeTypes <- c("MGW", "HelT")
shapeList <- getShape(tempFile, shapeType = shapeTypes)

write_rds(shapeList, paste0(outPrefix, ".shapeList.rds"))
# shapeList <- read_rds(paste0(outPrefix, ".shapeList.rds"))

names(shapeList)

# group anchors by participating in loop
ancLoop <- unique(unlist(anchors(gi[gi$loop == "Loop"], id = TRUE)))

# # convert in tidy tibble:
# shapeList <- map(shapeList, function(m) {
#   if (ncol(m) == SHAPE_WINDOW - 1 ) { m <- cbind(m, NA) } 
#   return(m)
# })

# summaryList <- map(shapeList, function(m) {
#   cmean <- colMeans(m, na.rm = TRUE)
#   cmedian <- colMedians(m, na.rm = TRUE)
#   sd <- colSds(m, na.rm = TRUE)
#   tibble(
#     pos = -50:50
#   )
#   return(m)
# })

#*******************************************************************************
# Plot shape around loop and non-loop motifs next to each other ----
#*******************************************************************************

pdf(paste0(outPrefix, ".shape_at_anchors_by_loop.pdf"), w = 6, h = 18)
par(mfrow = c(length(shapeTypes), 2))  # n rows and 2 columns
for (TYPE in shapeTypes) {
  
  print(TYPE)
  mat <- shapeList[[TYPE]]
  
  plotShape(mat[ancLoop, ], main = paste(TYPE, "at loop anchor"))
  plotShape(mat[-ancLoop, ], main = paste(TYPE, "at non-loop anchor"))
}
dev.off()

#*******************************************************************************
# Plot shape around loop and non-loop motifs together ----
#*******************************************************************************

pdf(paste0(outPrefix, ".shape_at_anchors_by_loop_together.pdf"), w = 12, h = 12)
par(mfrow = c(length(shapeTypes)/2, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,1,2,2) + 0.1)  # n rows and 2 columns

for (TYPE in shapeTypes) {
  
  print(TYPE)
  mat <- shapeList[[TYPE]]
  
  plotShape(
    shapeMatrix = mat[ancLoop, ], 
    background = mat[-ancLoop, ], 
    colLine = COL_LOOP[2], 
    colLineBg = COL_LOOP[1],
    main = TYPE)
  legend("topright", legend = c("Loop", "No loop"), 
         col = rev(COL_LOOP), lwd = 2 )
}
dev.off()

#*******************************************************************************
# Plot shape around loop and non-loop as heatmap  ----
#*******************************************************************************


pdf(paste0(outPrefix, ".heatShape_at_anchors_by_loop.pdf"), w = 6, h = 18)
par(mfrow = c(length(shapeTypes), 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,0,1,1) + 0.1)  # n rows and 2 columns

for (TYPE in shapeTypes) {
  
  print(TYPE)
  mat <- shapeList[[TYPE]]
  
  if (ncol(mat) == SHAPE_WINDOW) {
    mat <- cbind(mat, NA)
  }
  loopMat <- mat[ancLoop, ]
  nonLoopMat <- mat[-ancLoop, ]
  if( nrow(nonLoopMat) > 100){
    loopMat <- loopMat[sample.int(size = 100, n = nrow(loopMat)), ]
    nonLoopMat <- nonLoopMat[sample.int(size = 100, n = nrow(nonLoopMat)), ]
  }
  heatShape(loopMat, nBins = 20, main = paste(TYPE, "at loop anchor"))
  heatShape(nonLoopMat, nBins = 20, main = paste(TYPE, "at non-loop anchor"))
}

dev.off()





#*******************************************************************************
# Analyse correlation of shape features for motif pairs
#*******************************************************************************

# make a tidy DF
tidyDF <- mcols(gi) %>%
  as.data.frame() %>% 
  as.tibble() %>% 
  gather(starts_with("cor_"), key = name, value = cor) %>% 
  mutate(name = str_sub(name, 5)) 


p <- ggplot(tidyDF, aes(x = loop, y = cor)) +
  geom_violin(aes(fill = interaction(loop, name))) + 
  geom_boxplot(fill = "white", width = .2, outlier.shape = NA) +
  facet_grid(. ~ name) + 
  theme_bw() +
  scale_color_brewer(type = "qualitative", palette = "Set3") +
  theme(
    text = element_text(size = 20), 
    # axis.text.x=element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "DNA-shape\n Correlation", x = "")
ggsave(p, file = paste0(outPrefix, ".shape_cor.by_TF_and_loop.violin.pdf"), w = 7, h = 3.5)


p <- ggplot(tidyDF, aes(x = name, y = cor, col = loop)) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = COL_LOOP, guide_legend(title = "")) +
  theme_bw() + 
  theme(
    text = element_text(size=20),
    # axis.text.x=element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(y = "DNA-shape\n Correlation", x="")
ggsave(p, file = paste0(outPrefix, ".shape_cor.by_TF_and_loop.boxplot.pdf"), w = 7, h = 3.5)


