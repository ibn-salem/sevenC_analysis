#' fitts a logit model to training data and returns a tidy model fit.
#' 
#' @see broom::tidy()
#' 
#' @param splits an individual \code{rsplit} object.
#' @param design a formular object defining the prediction design and variables.
tidy_fitter <- function(splits, formula, ...) {
  
  # fit glm model
  mod <- glm(formula = formula, 
             data = analysis(splits), 
             family = binomial(),
             model = FALSE,  
             x = FALSE,
             ...)
  
  return(broom::tidy(mod))
}

#' Fitts a logit model to training data and returns a tidy model fit.
#' 
#' @see broom::tidy()
#' 
#' @param fold An indicator of current Fold.
#' @param tidyCV A tibble with columns Row, Data, and Fold. See \code{\link[rsample]{tidy.rsplit}}.
#' @param data a data.frame on which the cross-validation was executed.
#' @param design a formular object defining the prediction design and variables.
#' 
#' @value A tidy glm model. See \code{\link[broom]{tidy}}.
tidyer_fitter <- function(fold, formula, data, tidyCV, ...) {
  
  # extract row indices 
  rows <- tidyCV %>% 
    dplyr::filter(Data == "Analysis" & Fold == fold) %>% 
    pull("Row")
  
  # fit glm model only only the subset of rows in data
  mod <- glm(formula = formula, 
             data = dplyr::slice(data, rows), 
             family = binomial(),
             model = FALSE,  
             x = FALSE,
             ...)
  
  return(broom::tidy(mod))
}


#' Get assesment part of data based on cross-validation split using
#' \code{rsample}.
#'
#' @see \code\link[rsample]{assessment}}
#'
#' @param fold An indicator of current Fold.
#' @param tidyCV A tibble with columns Row, Data, and Fold. See
#'   \code{\link[rsample]{tidy.rsplit}}.
#' @param data a data.frame on which the cross-validation was executed.
#' 
#' @value The assessment subset of \code{data} according to fold.
tidy_assessment <- function(fold, data, tidyCV){
  
  # extract row indices 
  rows <- tidyCV %>% 
    dplyr::filter(Data == "Assessment" & Fold == fold) %>% 
    dplyr::pull("Row")
  
  # return subset of data
  return(dplyr::slice(data, rows))
}


#' Write interactions as long-range text file format.
#'
#'
#' The WashU EpiGenome Browser can visuallize chromaint interactions along other
#' genomic features. The format of the file format is described here:
#' http://wiki.wubrowse.org/Long-range
#'
#' @param gi
#' @param score_vec Numeric vector of the same length as \cod{gi}. It will be
#'   used as score in output fiel.
#' @param output_file Path to output file.
#'   
writeLongRangeFormat <- function(gi, score_vec, output_file){
  
  chr <- seqnames(regions(gi))
  s <- start(regions(gi))
  e <- end(regions(gi))
  
  # get indexes of anchor regions
  idx1 <- anchors(gi, type = "first", id = TRUE)
  idx2 <- anchors(gi, type = "second", id = TRUE)
  
  # gr1 <- anchors(gi, "first")
  # gr2 <- anchors(gi, "second")
  # # score_vec <- mcols(gi)[, score]
  
  # str1 <- str_c(seqnames(gr1), start(gr1), end(gr1), sep = ",")
  # str2 <- str_c(seqnames(gr2), start(gr2), end(gr2), sep = ",")
  str1 <- str_c(chr[idx1], s[idx1], e[idx1], sep = ",")
  str2 <- str_c(chr[idx2], s[idx2], e[idx2], sep = ",")

  outDF <- tibble(
    gr1 = str1, 
    gr2 = str2, 
    score = score_vec)
  
  write_tsv(outDF, path = output_file, col_names = FALSE)
}



#' Get contingency table as tidy tibble from.
#'
#' @param true Logical vector of true labels (goldstandard)
#' @param pred Logical vector of predicted lables
#' @return A tibble with the number of true postivies (TP), false postivives
#'   (FP), false negatives (FN), and true negatives (TN) in columns.
getContab <- function(true, pred){
  
  tab <- count(
    tibble(
      true = true,
      pred = pred
    ), 
    true, pred
  ) %>% 
    arrange(true, pred)
  
  tibble(
    TN = tab$n[[1]],
    FP = tab$n[[2]],
    FN = tab$n[[3]],
    TP = tab$n[[4]]
  )
}

# helper function
isLinked <- function(linkDF, len){
  # 1:len %in% linkDF$gr1
  is.element(seq(len), linkDF$gr1)
}


