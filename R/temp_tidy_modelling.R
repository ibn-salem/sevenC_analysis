

#' Get prediciton result on test subset
#' Source: https://github.com/topepo/rsample/blob/master/vignettes/Working_with_rsets.Rmd
#' 
#' @param splits an individual \code{rsplit} object.
#' 
holdout_results <- function(splits, formula, dependant_var = "loop", ...) {
  
  # Fit the model to the training data set
  mod <- glm(formula = formula, 
             data = analysis(splits), 
             family = binomial(),
             model = FALSE,  
             x = FALSE,
             ...)
  
  # `augment` will save the predictions with the holdout data set
  res <- broom::augment(
    mod, 
    newdata = assessment(splits), 
    type.predict = "response"
  ) %>% 
    as_tibble() %>% 
    select(label = one_of(dependant_var), pred = ".fitted")
  
  
  # Return the assessment data set with the additional columns
  return(res)
}


################################################################################
# debugging evalmod() with multiple models and datasets
library(precrec)

scores <- list(
  rnorm(20),
  rnorm(20),
  rnorm(20),
  rnorm(20)
)

lab1 <- sample(0:1, 20, replace = TRUE)
lab2 <- sample(0:1, 20, replace = TRUE)

labels <- list(
  lab1,
  lab2,
  lab1,
  lab2
)

crvs <- evalmod(
  scores = scores,
  labels = labels,
  modnames = c("A", "A", "B", "B"),
  dsids = c(1, 2, 1, 2)
)

autoplot(crvs)
crvsDF <- as.data.frame(crvs)

unique(crvsDF$modname)

#----------------
samps4 <- create_sim_samples(2, 10, 10, c("random", "excel"))

################################################################################

#######
library(plyr)
library(precrec)
require(ggplot2)

df <- data.frame(
  score1 = rnorm(50),
  score2 = rnorm(50),
  label = rep(c("pos", "neg"), c(25, 25)),
  fold = sample(cut(seq(1, 50), breaks = 5, labels = FALSE))
)

# Use plyr::daply to split df
s1 <- daply(df, .(fold), function(x)return(x$score1))
s2 <- daply(df, .(fold), function(x)return(x$score2))
df$label2 <- as.numeric(factor(df$label))
lb <- daply(df, .(fold), function(x)return(x$label2))

# Use precrec::join_scores and precrec::join_labels to format input data
scores <- join_scores(s1, s2, byrow = TRUE)
labs <- join_labels(lb, lb, byrow = TRUE)

# Calculate ROC and precision-recall
crvs <- evalmod(scores = scores, labels = labs, modnames =
                  rep(c("score1", "score2"), each=5), dsids = rep(1:5, times=2))

# get data.frame of curves
curveDF <- as.data.frame(crvs)
# check number of model names
table(curveDF$modname)

# plot using base plotting
plot(crvs)

# plot using ggplot2
autoplot(crvs)

sessionInfo()
#######

glm_coefs_rec <- function(splits, rec, ...) {
  
  # Estimate the parameters using the analysis data
  trained_rec <- prep(rec, training = analysis(splits), 
                      retain = TRUE,
                      verbose = FALSE)
  
  # Apply the transformations to the data set, save the 
  # predictor values, and convert to a matrix
  design_matrix <- juice(trained_rec, all_predictors())
  design_matrix <- as.matrix(design_matrix)
  
  # Get the outcome values and fit the model using `lm.fit`
  y <- juice(trained_rec, Attrition)
  # Convert a one column tible to a vector
  y <- getElement(y, "Attrition")
  
  mod <- glm.fit(x = design_matrix, y = y, family = binomial(), ...)
  as.data.frame(t(mod$coefficients))
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#' Simplify model by deliting havy and not used stuff.
#' Source: https://www.r-bloggers.com/trimming-the-fat-from-glm-models-in-r/
stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}



#-------------------------------------------------------------------------------



#' fitts a logit model to training data and returns prediction response vector
#' on test data
#' 
#' @param training_data data set with training data.
#' @param design a formular object defining the prediction design and variables.
fitter <- function(training_data, design) {
  glm(
    design, 
    family = binomial(link = 'logit'),
    data = training_data,
    model = FALSE, 
    x = FALSE
  ) 
}

#-------------------------------------------------------------------------------



require(biglm)
bm <- biglm::bigglm(design, data = as_tibble(cv$train[[1]]), family = binomial(link = 'logit'))
broom::augment(x=fm1, newdata = Data, type.predict = "response")


