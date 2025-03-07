#' Choose the prior of AMSSB using cross-validation.
#' 
#' This function does k-fold cross-validation for AMSSB and chooses the prior with the smallest cv error. Gaussian prior, Laplace prior and Cauchy prior are considered in this function.
#' 
#' @param X An n by G matrix of covariates.
#' @param Y The outcome vector with length n, or the n by 1 outcome matrix.
#' @param nfolds The total number of folds for cross-validation. Default is 5.
#' @param degree The degree of the spline basis. Default is 5.
#' @param knots The knots of the spline basis. Default is NULL.
#' @param Boundary.knots The boundary knots of the spline basis. Default is NULL.
#' @param spline_der_length The length of the spline second derivative. Default is 1000.
#' @param svd_threshold The threshold for the singular value decomposition. Default is 0.01.
#' @param loss Loss used as cross-validation error. Valid options are "L1" and "L2". Default is "L2".
#' @param verbose A logical value indicating whether to print the progress. Default is TRUE.
#' @param ... Other arguments used in AMSSB.
#' @return An object with S3 class "AMSSB" selected by cross-validation.
#' \tabular{ll}{
#'   \code{data} \tab The data information. \cr
#'  \tab \cr
#'  \code{beta} \tab The sparse estimation of coefficients. \cr
#' \tab \cr
#' \code{intercept} \tab The fitted value of intercept. \cr
#' \tab \cr
#' \code{sigma_noise} \tab The fitted value of noise level \eqn{\tilde{\sigma}}.\cr
#' \tab \cr
#' \code{mu} \tab The posterior mean of the coefficients. \cr
#' \tab \cr
#' \code{Sigma} \tab The posterior covariance matrix of the coefficients. \cr
#' \tab \cr
#' \code{gamma} \tab The posterior inclusion probability for each groups.\cr
#' \tab \cr
#' \code{selected_groups} \tab The positive group indexes selected by the algorithm.\cr
#' \tab \cr
#' \code{lambda} \tab The value of the hyperparameter \eqn{\lambda}. If em=FALSE, this value is the same as the input argument.\cr
#' \tab \cr
#' \code{prior} \tab The slab part of the coefficient prior.\cr
#' \tab \cr
#' \code{nu} \tab Degree of freedom when the slab is T distribution.\cr
#' }
#' @seealso \code{\link{AMSSB}}, \code{\link{predict.AMSSB}}
#' 
#' @examples
#' # load GVSSB library
#' library(GVSSB)
#' 
#' # generate covariates
#' n <- 200
#' G <- 300
#' X <- mvtnorm::rmvnorm(n, sigma=diag(G))
#'
#' # define additive model
#' f1 <- function(x) 5 * sin(x)
#' f2 <- function(x) 2 * (x^2 - 0.5)
#' f3 <- function(x) 2 * exp(x)
#' f4 <- function(x) 3 * x
#' 
#' # generate response vector
#' Y <- f1(X[,1]) + f2(X[,2]) + f3(X[,3]) + f4(X[,4]) + rnorm(n, 0, 1)
#' 
#' # fit AMSSB model with cross-validation
#' fit <- cv.AMSSB(X, Y, nfolds = 5, loss = 'L2')
#' 
#' 
#' @export cv.AMSSB
#' @export
cv.AMSSB = function(X, Y, nfolds = 5, degree = 5, knots = NULL, Boundary.knots = NULL, spline_der_length = 1000, svd_threshold = 0.01,  loss = 'L2', verbose = TRUE, ...){
  
  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  n = nrow(X)
  G = ncol(X)
  if(n != length(Y)) stop("The dimension of X and Y don't match!")
  if(!(loss %in% c('L2','L1'))) stop('Loss type incorrect. Only L1 loss and L2 loss are supported.')
  error = matrix(NA, nrow = nfolds, ncol = 3)
  Y = as.matrix(Y)

  # cross-validation
  for(fold in 1:nfolds){
    if(verbose) message(paste0('Fold ', fold, ' of ', nfolds, '...'))
    # separate training and testing data
        Xtrain = X[-(((fold- 1) * n / nfolds + 1):(fold * n / nfolds)), ]
        Ytrain = Y[-(((fold- 1) * n / nfolds+ 1):(fold * n / nfolds)), ]
        Xtest = X[(((fold - 1) * n / nfolds + 1):(fold * n / nfolds)), ]
        Ytest = Y[(((fold - 1) * n / nfolds + 1):(fold * n / nfolds)), ]

    # fit AMSSB model with priors
    fit.Gaussian = AMSSB(Xtrain, Ytrain, prior = 'Gaussian', degree = degree, knots = knots, Boundary.knots = Boundary.knots, spline_der_length = spline_der_length, svd_threshold = svd_threshold, info = F, ...)
    fit.Laplace = AMSSB(Xtrain, Ytrain, prior = 'Laplace', degree = degree, knots = knots, Boundary.knots = Boundary.knots, spline_der_length = spline_der_length, svd_threshold = svd_threshold, info = F, ...)
    fit.Cauchy = AMSSB(Xtrain, Ytrain, prior = 'T', nu=1, degree = degree, knots = knots, Boundary.knots = Boundary.knots, spline_der_length = spline_der_length, svd_threshold = svd_threshold, info = F, ...)

    # get the prediction
    predict.Gaussian = predict.AMSSB(fit.Gaussian, Xtest)
    predict.Laplace = predict.AMSSB(fit.Laplace, Xtest)
    predict.Cauchy = predict.AMSSB(fit.Cauchy, Xtest)

    # calculate the loss
    error[fold, 1] = loss_func(loss, Ytest, predict.Gaussian$Y_hat)
    error[fold, 2] = loss_func(loss, Ytest, predict.Laplace$Y_hat)
    error[fold, 3] = loss_func(loss, Ytest, predict.Cauchy$Y_hat)
  }

  # choose the best prior
  mean_error = colMeans(error)
  priors = c('Gaussian', 'Laplace', 'Cauchy')
  cv.prior = priors[which.min(mean_error)]
  if(verbose) message(paste0('The best prior is ', cv.prior, '. Fitting the model with full data...'))
  if(cv.prior == 'Cauchy') cv.prior = 'T'
  fit = AMSSB(X, Y, prior = cv.prior, degree = degree, knots = knots, Boundary.knots = Boundary.knots, spline_der_length = spline_der_length, svd_threshold = svd_threshold, info = F, ...)

  return(fit)

}
  