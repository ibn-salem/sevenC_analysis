

#' Extend interaction anchors with specific term for inner and outer site
#'
#' Note, this functions removes all metadata colums from the
#' \code{\link{GInteractions}} and its anchor regions.
#' 
#' 
#' @param gi \code{\link{GInteractions}} or \code{\link{InteractionSet}} object
#' @param inner,outer A scalar, non-negative, of the extension sizes for anchor
#'   ends outside and inside interaction loops.
#' @return gi \code{\link{GInteractions}} object with extened regions. The
#'   \code{\link{GInteractions}} object has two metadata columns, \code{dist},
#'   and \code{dist_log10} with the genomic distance between the midpointsof
#'   anchors and distance in log10 space, respectively.
extendAnchors <- function(gi, inner, outer){
  
  anchrosGR <- regions(gi)
  mcols(anchrosGR) <- NULL
  
  upAnc <- anchrosGR[anchors(gi, type = "first", id = TRUE)]
  downAnc <- anchrosGR[anchors(gi, type = "second", id = TRUE)]
  
  # extend upstream anchro range
  start(upAnc) <- start(upAnc) - outer
  end(upAnc) <- end(upAnc) + inner
  
  # extend downstream anchor range
  start(downAnc) <- start(downAnc) - inner
  end(downAnc) <- end(downAnc) + outer
  
  extGI <- InteractionSet::GInteractions(upAnc, downAnc)
  
  return(extGI)
}
  


#' Get genomic ranges spanned by (intracrhomosomal) interactions.
#'
#' Note, this method assumes are intra-chromosomal (e.g. only interactions
#' between regions from the same chromosome).
#'
#' @param gi \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object
#' @return A \code{\link[GenomicRanges]{GRanges}} object with ranges spanning
#'   the interactions in \code{gi}.
#'
#' @export
interactionRange <- function(gi){
  
  # assume only intra-chromosomal interactions
  stopifnot(all(InteractionSet::intrachr(gi)))
  
  # assume no redundant interactions
  stopifnot(length(gi) == length(swapAnchors(gi)))
  
  # swap anchors to make left anchor <= right anchor
  
  gi <- as(gi, "GInteractions")
  ancGR <- regions(gi)
  strand(ancGR) <- "*"
  regions(gi) <- ancGR
  gi <- swapAnchors(gi)

  # get al anchors as GRanges object
  ancGR <- regions(gi)
  
  # get chroms from first anchor
  upAncIdx <- anchors(gi, "first", id = TRUE)
  downAncIdx <- anchors(gi, "second", id = TRUE)
  
  # get chromosome, start of lef, and end of right anchor
  chr <- seqnames(ancGR)[upAncIdx]
  start <- start(ancGR)[upAncIdx]
  end <- end(ancGR)[downAncIdx]
  
  # assume start <- end
  stopifnot(all(start <= end))
  
  # create GenomicRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start, end),
    strand = "*",
    seqinfo = GenomeInfoDb::seqinfo(gi)
  )
  
  # add mcols
  S4Vectors::mcols(gr) <- S4Vectors::mcols(gi)
  
  return(gr)
}





#' Associate two sets of regions together by interactions and direct overlap
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param outer_maxgap A scalar, non-negative, maximal allowed distance outside
#'   interaction loops to be considered as overlap.
#' @param inner_maxgap A scalar, non-negative, maximal allowed distance inside
#'   interaction loops to be considered as overlap.
#' @param ... Other arguments passed to
#'   \code{\link[GenomicRanges]{findOverlaps}} and
#'   \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A tibble of integer indices indicating which elements of \code{gr1}
#'   link which elements of \code{gr1} by interaction \code{gi} (NA if direct
#'   overlap).
#'
#' @export
linkRegions <- function(gr1, gr2, gi, outer_maxgap = 0, inner_maxgap = 0, ...){

  # build new GInteraction object with extened anchor regions
  gi <- extendAnchors(gi, outer = outer_maxgap, inner = inner_maxgap)

  # get direct overlapping elements
  ol <- GenomicRanges::findOverlaps(gr1, gr2, ...)

  # link regions by interactions
  links <- InteractionSet::linkOverlaps(gi, gr1, gr2, ...)

  # take only gr1 and gr2 indices and rename to match direct overlap
  olDF <- ol %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    rename(gr1 = queryHits, gr2 = subjectHits) %>% 
    mutate(gi = NA)
  
  linksDF <- links %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    rename(gi = query, gr1 = subject1, gr2 = subject2)
    
  hits <- bind_rows(olDF, linksDF)
  
  return(hits)
}


#' Associate two sets of regions together when they are within a loop.
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param ... Other arguments passed to \code{\link[GenomicRanges]{findOverlaps}}
#'   and \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A dataframe of integer indices indicating which elements of
#'   \code{gr1} link which elements of \code{gr1}.
#'
#' @export
linkRegionsInLoops <- function(gr1, gr2, gi, ...){

  # get range of all looping interactions
  loopGR <- interactionRange(gi)

  gr1Hits <- GenomicRanges::findOverlaps(gr1, loopGR) %>% 
    data.frame() %>% as.tibble() %>% 
    rename(gr1 = queryHits)

  gr2Hits <- GenomicRanges::findOverlaps(gr2, loopGR) %>% 
    data.frame() %>% as.tibble() %>% 
    rename(gr2 = queryHits)
  
  # gr2Hits <- data.frame(GenomicRanges::findOverlaps(gr2, loopGR))
  # gr2Hits <- dplyr::rename(gr2Hits, gr2 = queryHits)

  # link reginos by
  linksDF <- left_join(gr1Hits, gr2Hits, by = "subjectHits") %>% 
    select(gi = subjectHits, gr1, gr2)

  # get direct overlapping elements
  olDF <- GenomicRanges::findOverlaps(gr1, gr2, ...) %>% 
    as.data.frame() %>% as.tibble() %>% 
    mutate(gi = NA) %>% 
    select(gi, gr1 = queryHits, gr2 = subjectHits)
  
  # combine direct overlaps and links via interactions
  hits <- bind_rows(olDF, linksDF)

  return(hits)
}
