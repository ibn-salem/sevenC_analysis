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


#' Predict outcome using logistic regresstion.
#' 
#' Source: 
#' https://stackoverflow.com/questions/25695565/predict-with-arbitrary-coefficients-in-r
#' 
#' @param data  A data.frame with predictor variables
#' @param formula A modelling formula
#' @param betas A vector with parameter estimates for predictor variables
#'   
#' @return a data.frame like \cod{df} with additional columns with the 
#'   predictions.
#'   
#' @importFrom boot inv.logit
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate select
pred_logit <- function(data, formula, betas){
  
  # save option and set na.action to pass
  op <- options(na.action = 'na.pass')
  on.exit(options(op))
  
  # build model matrix
  mat <- modelr::model_matrix(data, formula)
  
  pred <- as.numeric(betas %*% t(mat)) %>%
    boot::inv.logit()
  
  return(pred)
  
}

# data <- assessment(dfCV$splits[[1]])
# formula <- dfCV$design[[1]]
# betas <- dfCV$tidy_model[[1]]$estimate
# pred <- pred_logit(data, formula, betas)

test_pred_logit <- function(){
  
  # testing
  formula <- am ~ wt + cyl + hp 
  mod <- glm(formula, mtcars, family = binomial())
  betas <- coef(mod)
  
  pred <- pred_logit(mtcars, betas, formula = formula)
  pred_glm <- predict.glm(mod, newdata = mtcars, type = "response")
  names(pred_glm) <- NULL
  
  # tibble(
  #   pred_glm = pred_glm, 
  #   pred = pred
  # )
  stopifnot(all.equal(pred, pred_glm))
}

test_pred_logit()


