################################################################################
# Analysis Motifs in true loop anchors
################################################################################

require(chromloop)    # for interaction parsers
require(InteractionSet) # for GI objects
require(tidyverse)    # for tidy data
require(stringr)      # for string functions
require(rtracklayer)  # to import() BED files
require(stringr)    # for string function and regular expressions
require(TxDb.Hsapiens.UCSC.hg19.knownGene)  # for seqinfo object
require(BSgenome.Hsapiens.UCSC.hg19)

# Parameters ===================================================================

outPrefix <- file.path("results", "motif_enrichment")
dir.create(outPrefix, showWarnings = FALSE)


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

# Parse loops ==================================================================
seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
genome <- BSgenome.Hsapiens.UCSC.hg19

trueLoopsRao <- chromloop::parseLoopsRao(
  LoopRao2014_GM12878_File, seqinfo = seqInfo)

trueLoopsTang2015 <- do.call(
  "c",
  lapply(LoopTang2015_GM12878_Files, 
         chromloop::parseLoopsTang2015, 
         seqinfo = seqInfo))

trueLoopsTangCTCF <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[1]], seqinfo = seqInfo)
trueLoopsTangPolII <- chromloop::parseLoopsTang2015(LoopTang2015_GM12878_Files[[2]], seqinfo = seqInfo)

trueCaptureHiC <- do.call(
  "c",
  lapply(CaptureHiC_Files,
         chromloop::parseCaptureHiC, seqinfo = seqInfo)
)


allGI <- list(
  HiC = trueLoopsRao,
  ChIAPet_CTCF = trueLoopsTangCTCF,
  ChIAPet_PolII = trueLoopsTangPolII,
  CaptureHiC = trueCaptureHiC
)

# Get anchors ==================================================================

ancGRL <- GRangesList(lapply(allGI, function(gi){
  sort(
    trim(
      # GenomicRanges::reduce(
        InteractionSet::regions(gi)
      # )
    )
  )
}))

ancGR <- sort(unlist(ancGRL))

# Get sequences and output as fasta file =======================================

ancSeq <- getSeq(genome, ancGR)

writeXStringSet(
    ancSeq, 
    filepath = file.path(outPrefix, "all_anchors.raw.fasta"), 
    format = "fasta")

for (name in names(ancGRL)) {
  
  message("INFO: working on ", name)
  
  gr <- ancGRL[[name]]
  
  seq <- getSeq(genome, gr)
  names(seq) <- str_c(name, "_", seq_along(seq))
  
  writeXStringSet(
    seq, 
    filepath = file.path(outPrefix, paste0(name, "_anchors.raw.fasta")), 
    format = "fasta")
}







