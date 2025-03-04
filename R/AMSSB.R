#' Fit an additive model with discrete spike and slab prior in variational Bayes.
#' 
#' This function gives a sparse estimation of the coefficients in additive model using spike and slab variational Bayes. To make prediction, please use \code{\link{predict.AMSSB}}.
#' 
#' @import splines2
#' 
#' @param X An n by G matrix of covariates.
#' @param Y The outcome vector with length n, or the n by 1 outcome matrix.
#' @param degree The degree of the spline basis. Default is 5.
#' @param knots The knots of the spline basis. Default is NULL, which means that the knots are equally spaced.
#' @param Boundary.knots The boundary knots of the spline basis. Default is NULL, which means that the boundary knots are the minimum and maximum of the covariates.
#' @param spline_der_length The total length of intervals for the spline second derivative. Default is 1000.
#' @param svd_threshold The threshold for the singular value decomposition of the kernel matrix. Any singular value below this threshold will be set to this threshold. Default is 0.01.
#' @param ... Other arguments used in GVSSB.
#' @return An object with S3 class "AMSSB".
#' \tabular{ll}{
#'    \code{data} \tab The data information. \cr
#'    \tab \cr
#'    \code{beta} \tab The sparse estimation of coefficients. \cr
#'    \tab \cr
#'    \code{intercept} \tab The fitted value of intercept. \cr
#'    \tab \cr
#'    \code{sigma_noise} \tab The fitted value of noise level \eqn{\tilde{\sigma}}.\cr
#'    \tab \cr
#'    \code{mu} \tab The posterior mean of the coefficients. \cr
#'    \tab \cr
#'    \code{Sigma} \tab The posterior covariance matrix of the coefficients. \cr
#'    \tab \cr
#'    \code{gamma} \tab The posterior inclusion probability for each groups.\cr
#'    \tab \cr
#'    \code{selected_groups} \tab The positive group indexes selected by the algorithm.\cr
#'    \tab \cr
#'    \code{lambda} \tab The value of the hyperparameter \eqn{\lambda}. If em=FALSE, this value is the same as the input argument.\cr
#'    \tab \cr
#'    \code{prior} \tab The slab part of the coefficient prior.\cr
#'    \tab \cr
#'    \code{nu} \tab Degree of freedom when the slab is T distribution.\cr
#' }
#' @seealso \code{\link{GVSSB}}, \code{\link{predict.AMSSB}}, \code{\link{cv.GVSSB}}
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
#' @export 
AMSSB = function(X, Y, degree = 5, knots = NULL, Boundary.knots = NULL, spline_der_length = 1000, svd_threshold = 0.01, ...){
  
  # consistency check
  if(!is.matrix(X)) stop("X should be a matrix!")
  n = nrow(X)
  G = ncol(X)
  p = G * degree
  groups = rep(1:G, each = degree)
  if(n != length(Y)) stop("The dimension of X and Y don't match!")
  if(!is.null(knots)){
    if(!is.matrix(knots)) stop("knots should be a matrix!")
    if(nrow(knots) != degree-3) stop("The number of knots should be degree-3!")
    if(ncol(knots) != G) stop("The number of knots should be the same as the number of covariates!")
  }
  if(!is.null(Boundary.knots)){
    if(!is.matrix(Boundary.knots)) stop("Boundary.knots should be a matrix!")
    if(nrow(Boundary.knots) != 2) stop("The number of Boundary.knots should be 2!")
    if(ncol(Boundary.knots) != G) stop("The number of Boundary.knots should be the same as the number of covariates!")
  }

  # standardize the covariates
  X.min = apply(X, 2, min)
  X.max = apply(X, 2, max)
  X.ran = X.max - X.min
  X.min.rep = matrix(rep(X.min,n),nrow=n,byrow=T)
  X.ran.rep = matrix(rep(X.ran,n),nrow=n,byrow=T)
  X.norm = (X-X.min.rep)/X.ran.rep

  # get spline expansion for each covariate
  Z = matrix(0,n,p)
  nknots_bs = matrix(0,degree-3,G)
  Boundary.knots_bs = matrix(0,2,G)
  for(j in 1:G){
    index = (j-1)*degree + c(1:degree)
    bsMat = bSpline(X.norm[,j],df=degree,knots=knots[,j],Boundary.knots=Boundary.knots[,j])
    Z[,index] = bsMat
    nknots_bs[,j] = attr(bsMat,'knots')
    Boundary.knots_bs[,j] = attr(bsMat,'Boundary.knots')
  }

  # get the inverse of the kernel matrix (spline second derivative)
  d_spline = matrix(0,p, spline_der_length)
  for(j in 1:G){
      index = (j-1)*degree + c(1:degree)
      x = seq(Boundary.knots_bs[1,j],Boundary.knots_bs[2,j],length.out=spline_der_length)
      d_spline[index,] = t(bSpline(x, derivs = 2, df =degree, knots = nknots_bs[,j], Boundary.knots = Boundary.knots_bs[,j]))
  }
  d_spline_kernel  = matrix(0,p,p)
  delta = (Boundary.knots_bs[2,]-Boundary.knots_bs[1,])/spline_der_length
  for(i in 1:G){
    index = (i-1)*degree + c(1:degree)
    d_spline_kernel[index,index] = d_spline[index,]%*%t(d_spline[index,])*delta[i]
  }
  d_spline_kernel_inv = matrix(0,p,p)
  for(i in 1:G){
    index = (i-1)*degree + c(1:degree)
    svd_kernel =  svd(d_spline_kernel[index,index])
    u_kernel = svd_kernel$u
    d_kernel = svd_kernel$d
    d_kernel[d_kernel < svd_threshold] = svd_threshold
    d_spline_kernel_inv[index,index] = solve(u_kernel%*%diag(d_kernel)%*%t(u_kernel))
  }

  # fit the model
  model = GVSSB(Z, Y, groups, standardize = F, Omega = d_spline_kernel_inv, ...)
  data_model = list(X.min = X.min, X.max = X.max, degree = degree, knots = nknots_bs, Boundary.knots = Boundary.knots_bs, Omega = d_spline_kernel_inv)
  model$data = data_model
  class(model) = 'AMSSB'
  return(model)
}