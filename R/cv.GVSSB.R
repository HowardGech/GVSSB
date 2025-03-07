#' Choose the prior of GVSSB using cross-validation.
#'
#' This function does k-fold cross-validation for GVSSB and chooses the prior with the smallest cv error. Gaussian prior, Laplace prior and Cauchy prior are considered in this function.
#'
#' @param X An n by p matrix of covariates.
#' @param Y The outcome vector with length n, or the n by 1 outcome matrix.
#' @param groups The group indicator vector with length p.
#' @param nfolds The total number of folds for cross-validation. Default is 5.
#' @param loss Loss used as cross-validation error. Valid options are "L1" and "L2". Default is "L2".
#' @param verbose A logical value indicating whether to print the progress. Default is TRUE.
#' @param ... Other arguments used in GVSSB.
#' @return An object with S3 class "GVSSB" selected by cross-validation.
#' \tabular{ll}{
#'    \code{beta} \tab The sparse estimation of coefficients. \cr
#'    \tab \cr
#'    \code{intercept} \tab The fitted value of intercept. \cr
#'    \tab \cr
#'    \code{sigma} \tab The fitted value of noise level \eqn{\tilde{\sigma}}.\cr
#'    \tab \cr
#'    \code{gamma} \tab The posterior inclusion probability for each groups.\cr
#'    \tab \cr
#'    \code{lambda} \tab The value of the hyperparameter \eqn{\lambda}. If em=FALSE, this value is the same as the input argument.\cr
#'    \tab \cr
#'    \code{groups_index} \tab The group index with length p.\cr
#'    \tab \cr
#'    \code{prior} \tab The slab part of the coefficient prior.\cr
#'    \tab \cr
#'    \code{nu} \tab Degree of freedom when the slab is T distribution.\cr
#' }
#' @seealso \code{\link{GVSSB}}, \code{\link{predict.GVSSB}} 
#' @examples
#' # load GVSSB library
#' library(GVSSB)
#' 
#' # generate covariates
#' n <- 200
#' G <- 200
#' p_i <- 5
#' p <- G * p_i
#' X <- mvtnorm::rmvnorm(n, sigma=diag(p))
#' 
#' # generate coefficients
#' k <- 5
#' beta <- rep(0,p)
#' nonzero_group <- sample(1:G, k)
#' for(index in nonzero_group){
#'     beta[p_i * (index - 1) + 1:p_i] <- runif(p_i, -1, 1)
#' }
#' 
#' # define group index
#' groups <- rep(1:G, each = p_i)
#' 
#' # generate response vector
#' Y <- X %*% beta + rnorm(n, 0, sd = 1)
#' 
#' # fit GVSSB model with cross-validation
#' fit <- cv.GVSSB(X, Y, groups, nfolds = 5, loss = 'L2')
#' 
#' 
#' @export cv.GVSSB
#' @export
cv.GVSSB = function(X, Y, groups, nfolds = 5, loss = c('L2','L1'), verbose = TRUE, ...){
  loss = match.arg(loss)
  Y = as.matrix(Y)
  n = nrow(Y)

  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  if(!(loss %in% c('L2','L1'))) stop('Loss type incorrect. Only L1 loss and L2 loss are supported.')
  if(nrow(X) != n) stop("The dimension of X and Y don't match!")

  error = matrix(NA, nrow = 5, ncol = 3)
  
  # cross-validation
  for(fold in 1:nfolds){
    if(verbose) message(paste0('Fold ', fold, ' of ', nfolds, '...'))
    # separate training data and testing data
    Xtrain = X[-(((fold- 1) * n / nfolds + 1):(fold * n / nfolds)), ]
    Ytrain = Y[-(((fold- 1) * n / nfolds+ 1):(fold * n / nfolds)), ]
    Xtest = X[(((fold - 1) * n / nfolds + 1):(fold * n / nfolds)), ]
    Ytest = Y[(((fold - 1) * n / nfolds + 1):(fold * n / nfolds)), ]
    
    # fit the GVSSB model with three priors
    fit.Gaussian = GVSSB(Xtrain, Ytrain, groups, prior = "Gaussian", info = F, ...)
    fit.Laplace = GVSSB(Xtrain, Ytrain, groups, prior = "Laplace", info = F, ...)
    fit.Cauchy = GVSSB(Xtrain, Ytrain, groups, prior = "T", info = F, ...)
    
    # get the predicted response vector
    predict.Gaussian = predict(fit.Gaussian, Xtest)
    predict.Laplace = predict(fit.Laplace, Xtest)
    predict.Cauchy = predict(fit.Cauchy, Xtest)
    
    # compute the loss
    error[fold, 1] = loss_func(loss, Ytest, predict.Gaussian$Y_hat)
    error[fold, 2] = loss_func(loss, Ytest, predict.Laplace$Y_hat)
    error[fold, 3] = loss_func(loss, Ytest, predict.Cauchy$Y_hat)
  }
  
  # choose the prior with the smallest cv error
  mean_error = colMeans(error)
  priors = c('Gaussian','Laplace','T')
  cv_prior = priors[which.min(mean_error)]
  if(verbose) message(paste0('The best prior is ', cv_prior, '. Fitting the model with full data...'))
  if(cv.prior == 'Cauchy') cv.prior = 'T'
  fit = GVSSB(X, Y, groups, prior = cv_prior, info = F, ...)

return(fit)
}