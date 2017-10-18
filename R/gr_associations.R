

#' Extend interaction anchors with specific term for inner and outer site
#'
#' @param gi \code{\link{GInteractions}} or \code{\link{InteractionSet}} object
#' @param inner,outer A scalar, non-negative, of the extension sizes for anchor
#'   ends outside and inside interaction loops.
#' @return gi \code{\link{GInteractions}} object with extened regions. The
#'   \code{\link{GInteractions}} object has two metadata columns, \code{dist},
#'   and \code{dist_log10} with the genomic distance between the midpointsof
#'   anchors and distance in log10 space, respectively.
extendAnchors <- function(gi, inner, outer){
  
  # turn gi into GIntreactions
  gi <- methods::as(gi, "GInteractions")
  GenomicRanges::strand(InteractionSet::regions(gi)) <- "*"
  gi <- InteractionSet::swapAnchors(gi, mode = "order")
  
  # extend ranges of anchors
  anc1 <- InteractionSet::anchors(gi, "first")
  GenomicRanges::start(anc1) = GenomicRanges::start(anc1) - outer
  # extend end coordinate of fist anchor but only until start of second
  GenomicRanges::end(anc1) = pmin(
    GenomicRanges::end(anc1) + inner,
    pmax(GenomicRanges::start(InteractionSet::anchors(gi, "second")) - 1, GenomicRanges::end(anc1))
  )
  
  anc2 <- InteractionSet::anchors(gi, "second")
  # extend start coordinate of second anchor but only until end of first
  GenomicRanges::start(anc2) = pmax(
    GenomicRanges::start(anc2) - inner,
    pmin(GenomicRanges::end(InteractionSet::anchors(gi, "first")) + 1, GenomicRanges::start(anc2))
  )
  GenomicRanges::end(anc2) = GenomicRanges::end(anc2) + outer
  
  extendedGI <- InteractionSet::GInteractions(anc1, anc2)
  S4Vectors::mcols(extendedGI) <- S4Vectors::mcols(gi)
  
  return(extendedGI)
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
  
  # make first anchor <= second anchor
  gi <- methods::as(gi, "GInteractions")
  GenomicRanges::strand(InteractionSet::regions(gi)) <- "*"
  gi <- InteractionSet::swapAnchors(gi, mode = "order")
  
  # get chroms from first anchor
  chr <- GenomeInfoDb::seqnames(InteractionSet::anchors(gi, "first"))
  
  # take all coordinates to get smalles as start and largest as end
  coords <- list(
    GenomicRanges::start(InteractionSet::anchors(gi, "first")),
    GenomicRanges::end(InteractionSet::anchors(gi, "first")),
    GenomicRanges::start(InteractionSet::anchors(gi, "second")),
    GenomicRanges::end(InteractionSet::anchors(gi, "second"))
  )
  start <- do.call(pmin, coords)
  end <- do.call(pmax, coords)
  
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





#' Associate two sets of regions together by interactions with ajusted anchors
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param outer_maxgap A scalar, non-negative, maxial allowed distance outside
#'   interaction loops to be considered as overlap.
#' @param inner_maxgap A scalar, non-negative, maxial allowed distance inside
#'   interaction loops to be considered as overlap.
#' @param ... Other arguments passed to
#'   \code{\link[GenomicRanges]{findOverlaps}} and
#'   \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A dataframe of integer indices indicating which elements of
#'   \code{gr1} link which elements of \code{gr1}.
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
  olDF <- as.data.frame(ol)
  names(olDF) <- c("gr1", "gr2")
  olDF[, "gi"] <- NA

  linksDF <- data.frame(links)
  names(linksDF) <- c("gi", "gr1", "gr2")

  # combine direct overlaps and links via interactions
  hits <- rbind(
    olDF,
    linksDF
  )

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

  # get direct overlapping elements
  ol <- GenomicRanges::findOverlaps(gr1, gr2, ...)

  # get range of all looping interactions
  loopGR <- interactionRange(gi)

  gr1Hits <- data.frame(GenomicRanges::findOverlaps(gr1, loopGR))
  gr1Hits <- dplyr::rename(gr1Hits, gr1 = queryHits)

  gr2Hits <- data.frame(GenomicRanges::findOverlaps(gr2, loopGR))
  gr2Hits <- dplyr::rename(gr2Hits, gr2 = queryHits)

  # link reginos by
  linksDF <- dplyr::left_join(gr1Hits, gr2Hits, by = "subjectHits")
  linksDF <- dplyr::select(linksDF, gi = subjectHits, gr1, gr2)


  # take only gr1 and gr2 indices and rename to match direct overlap
  olDF <- data.frame(ol)
  names(olDF) <- c("gr1", "gr2")
  olDF[, "gi"] <- NA

  # combine direct overlaps and links via interactions
  hits <- rbind(
    olDF,
    linksDF
  )

  return(hits)
}
