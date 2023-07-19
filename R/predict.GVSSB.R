#' Gives a prediction of a linear regression problem using a fitted GVSSB model.
#'
#' This function takes in a fitted GVSSB model and a covariate matrix, and gives a prediction of the response vector.
#'
#' @param object A fitted GVSSB object.
#' @param X An n by p matrix of covariates.
#' @param ... Other arguments used in GVSSB.
#' @return An object with S3 class "predicted_GVSSB".
#' \tabular{ll}{
#'    \code{Y_hat} \tab The estimated value of the response vector. \cr
#' }
#' @seealso \code{\link{GVSSB}}, \code{\link{cv.GVSSB}} 
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
#' k <- 10
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
#' Y <- X %*% beta + rnorm(n, 0, 1)
#' 
#' # fit GVSSB model
#' fit.Gaussian <- GVSSB(X, Y, groups, prior = 'Gaussian')
#' fit.Laplace <- GVSSB(X, Y, groups, prior = 'Laplace')
#' fit.Cauchy <- GVSSB(X, Y, groups, prior = 'T', nu = 1)
#' 
#' # 
#' @export predict.GVSSB
#' @export
predict.GVSSB = function(object, X, threshold, ...){
  
  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  n = nrow(X)
  p = ncol(X)
  if(p != length(object$groups_index)) stop("The dimension of X doesn't match!")
  
  # get the coefficients estimated by the GVSSB model
  beta = object$beta
  
  intercept = object$intercept
  
  # compute the prediction of the response vector
  Y_hat = X %*% beta + intercept
  fit = list(Y_hat = Y_hat)
  class(fit) = 'predict_GVSSB'
  return(fit)
}