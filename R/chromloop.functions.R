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
    filter(Data == "Analysis" & Fold == fold) %>% 
    pull("Row")
  
  # fit glm model only only the subset of rows in data
  mod <- glm(formula = formula, 
             data = slice(data, rows), 
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
    filter(Data == "Assessment" & Fold == fold) %>% 
    pull("Row")
  
  # return subset of data
  return(slice(data, rows))
}



