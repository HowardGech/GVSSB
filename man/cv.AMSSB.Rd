% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.AMSSB.R
\name{cv.AMSSB}
\alias{cv.AMSSB}
\title{Choose the prior of AMSSB using cross-validation.}
\usage{
cv.AMSSB(
  X,
  Y,
  nfolds = 5,
  degree = 5,
  knots = NULL,
  Boundary.knots = NULL,
  spline_der_length = 1000,
  svd_threshold = 0.01,
  loss = "L2",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{X}{An n by G matrix of covariates.}

\item{Y}{The outcome vector with length n, or the n by 1 outcome matrix.}

\item{nfolds}{The total number of folds for cross-validation. Default is 5.}

\item{degree}{The degree of the spline basis. Default is 5.}

\item{knots}{The knots of the spline basis. Default is NULL.}

\item{Boundary.knots}{The boundary knots of the spline basis. Default is NULL.}

\item{spline_der_length}{The length of the spline second derivative. Default is 1000.}

\item{svd_threshold}{The threshold for the singular value decomposition. Default is 0.01.}

\item{loss}{Loss used as cross-validation error. Valid options are "L1" and "L2". Default is "L2".}

\item{verbose}{A logical value indicating whether to print the progress. Default is TRUE.}

\item{...}{Other arguments used in AMSSB.}
}
\value{
An object with S3 class "AMSSB" selected by cross-validation.
\tabular{ll}{
  \code{data} \tab The data information. \cr
 \tab \cr
 \code{beta} \tab The sparse estimation of coefficients. \cr
\tab \cr
\code{intercept} \tab The fitted value of intercept. \cr
\tab \cr
\code{sigma_noise} \tab The fitted value of noise level \eqn{\tilde{\sigma}}.\cr
\tab \cr
\code{mu} \tab The posterior mean of the coefficients. \cr
\tab \cr
\code{Sigma} \tab The posterior covariance matrix of the coefficients. \cr
\tab \cr
\code{gamma} \tab The posterior inclusion probability for each groups.\cr
\tab \cr
\code{selected_groups} \tab The positive group indexes selected by the algorithm.\cr
\tab \cr
\code{lambda} \tab The value of the hyperparameter \eqn{\lambda}. If em=FALSE, this value is the same as the input argument.\cr
\tab \cr
\code{prior} \tab The slab part of the coefficient prior.\cr
\tab \cr
\code{nu} \tab Degree of freedom when the slab is T distribution.\cr
}
}
\description{
This function does k-fold cross-validation for AMSSB and chooses the prior with the smallest cv error. Gaussian prior, Laplace prior and Cauchy prior are considered in this function.
}
\examples{
# load GVSSB library
library(GVSSB)

# generate covariates
n <- 200
G <- 300
X <- mvtnorm::rmvnorm(n, sigma=diag(G))

# define additive model
f1 <- function(x) 5 * sin(x)
f2 <- function(x) 2 * (x^2 - 0.5)
f3 <- function(x) 2 * exp(x)
f4 <- function(x) 3 * x

# generate response vector
Y <- f1(X[,1]) + f2(X[,2]) + f3(X[,3]) + f4(X[,4]) + rnorm(n, 0, 1)

# fit AMSSB model with cross-validation
fit <- cv.AMSSB(X, Y, nfolds = 5, loss = 'L2')


}
\seealso{
\code{\link{AMSSB}}, \code{\link{predict.AMSSB}}
}
