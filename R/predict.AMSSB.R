#' Gives the prediction of an additive model problem using a fitted AMSSB model.
#' 
#' This function takes in a fitted AMSSB model and a covariate matrix, and gives a prediction of the response vector.
#' 
#' @import splines2
#' 
#' @param object A fitted AMSSB object.
#' @param X An n by G matrix of covariates.
#' @param truncate A logical value indicating whether to truncate the covariates. Default is TRUE.
#' @param ... Other arguments used in AMSSB.
#' @return An object with S3 class "predicted_AMSSB".
#' \tabular{ll}{
#'    \code{X} \tab The covariate matrix. \cr
#'    \tab \cr
#'    \code{Y_hat} \tab The estimated value of the response vector. \cr
#'    \tab \cr
#'    \code{model} \tab The fitted AMSSB model. \cr
#' }
#' @seealso \code{\link{AMSSB}}, \code{\link{cv.AMSSB}}
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
#' # fit AMSSB model
#' model.Laplace <- AMSSB(X, Y, prior = 'Laplace')
#' 
#' # generate new data
#' n_test <- 50
#' X_test <- mvtnorm::rmvnorm(n_test, sigma=diag(G))
#' 
#' # get the prediction
#' predict.Laplace <- predict.AMSSB(model.Laplace, X_test)
#' 
#' @export predict.AMSSB
#' @export
predict.AMSSB = function(object, X, truncate = TRUE, ...){

  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  n = nrow(X)
  G = ncol(X)
  X.min = object$data$X.min
  X.max = object$data$X.max
  if(G != length(X.min)) stop("The dimension of X and the fitted model don't match!")

  # standardize the covariates
  X.ran = X.max - X.min
  X.min.rep = matrix(rep(X.min, each = n), nrow = n, byrow = TRUE)
  X.ran.rep = matrix(rep(X.ran, each = n), nrow = n, byrow = TRUE)
  X_trans = (X - X.min.rep) / X.ran.rep

  if(truncate){
    X_trans = pmax(X_trans, 0)
    X_trans = pmin(X_trans, 1)
  }

  # get the spline expansion for each covariate
  knots = object$data$knots
  Boundary.knots = object$data$Boundary.knots
  degree = object$data$degree
  Z = matrix(0, n, G * degree)
  for(g in 1:G){
    index = (g - 1) * degree + c(1:degree)
    Z[, index] = bSpline(X_trans[,g],df=degree,knots=knots[,g],Boundary.knots=Boundary.knots[,g])
  }


  # get the coefficients estimated by the GVSSB model
  beta = object$beta

  intercept = object$intercept

  # compute the prediction of the response vector
  Y_hat = Z %*% beta + intercept
  fit = list(X = X, Y_hat = Y_hat, model = object)
  class(fit) = 'predicted_AMSSB'
  return(fit)
}
