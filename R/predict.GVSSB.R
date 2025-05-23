#' Gives the prediction of a linear regression problem using a fitted GVSSB model.
#'
#' This function takes in a fitted GVSSB model and a covariate matrix, and gives a prediction of the response vector.
#'
#' @param object A fitted GVSSB object.
#' @param X An n by p matrix of covariates.
#' @param ... Other arguments used in GVSSB.
#' @return An object with S3 class "predicted_GVSSB".
#' \tabular{ll}{
#'    \code{X} \tab The covariate matrix. \cr
#'    \tab \cr
#'    \code{Y_hat} \tab The estimated value of the response vector. \cr
#'    \tab \cr
#'    \code{model} \tab The fitted GVSSB model. \cr
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
#' Y <- X %*% beta + rnorm(n, 0, 1)
#'
#' # fit GVSSB model
#' fit.Gaussian <- GVSSB(X, Y, groups, prior = 'Gaussian')
#' fit.Laplace <- GVSSB(X, Y, groups, prior = 'Laplace')
#' fit.Cauchy <- GVSSB(X, Y, groups, prior = 'T', nu = 1)
#'
#' # generate new data
#' n_test <- 50
#' X_test <- mvtnorm::rmvnorm(n_test, sigma=diag(p))
#' 
#' # get the prediction
#' predict.Gaussian <- predict.GVSSB(fit.Gaussian, X_test)
#' predict.Laplace <- predict.GVSSB(fit.Laplace, X_test)
#' predict.Cauchy <- predict.GVSSB(fit.Cauchy, X_test)
#'
#' @export predict.GVSSB
#' @export
predict.GVSSB = function(object, X, ...){

  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  n = nrow(X)
  p = ncol(X)
  if(p != length(object$data$groups)) stop("The dimension of X doesn't match!")

  # get the coefficients estimated by the GVSSB model
  beta = object$beta

  intercept = object$intercept

  # compute the prediction of the response vector
  Y_hat = X %*% beta + intercept
  fit = list(X = X, Y_hat = Y_hat, model = object)
  class(fit) = 'predicted_GVSSB'
  return(fit)
}
