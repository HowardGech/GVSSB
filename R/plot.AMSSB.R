#' Gives the spline fitting plot of an additive model by AMSSB
#' 
#' This function takes in a fitted AMSSB model and the group index to be plotted, and generates a spline fitting plot based on the fitted values.
#' 
#' @import splines2
#' 
#' @param object A fitted AMSSB object.
#' @param index The index of the covariate to be plotted.
#' @param subtract_mean A logical value indicating whether to subtract the mean from the fitted values. Default is TRUE.
#' @param length_out The number of points to be used for the spline fitting.
#' @param ... Additional arguments to be passed to the plot function.
#'
#' @return A plot of the spline fitting for the specified group.
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
#' # plot spline fitting for group 1
#' plot.AMSSB(model.Laplace, index = 1)
#' 
#' @export plot.AMSSB
#' @export
plot.AMSSB = function(object, index = NULL, subtract_mean = TRUE, length_out = 500, ...) {
  # Check if the object is of class AMSSB
  if (!inherits(object, "AMSSB")) {
    stop("The object must be of class 'AMSSB'.")
    }

    selected_groups = object$selected_groups
    if(length(selected_groups) == 0) {
      stop("No groups selected by the model.")
    }
    if (is.null(index)) {
      index = selected_groups[1]
    } else {
      if (!index %in% selected_groups) {
        stop("The specified group index is not in the selected groups.")
      }
    }
    X.min = object$data$X.min
    X.max = object$data$X.max
    X.seq = seq(X.min[index], X.max[index], length.out = length_out)
    degree = object$data$degree
    data.seq = seq(0, 1, length.out = length_out)
    G = length(X.min)
    X_arr = matrix(0, nrow = length_out, ncol = G*degree)
    knots = object$data$knots
    Boundary.knots = object$data$Boundary.knots
    group_index = (index - 1) * degree + 1:degree
    X_arr[, group_index] = bSpline(data.seq,df=degree,knots=knots[,index],Boundary.knots=Boundary.knots[,index])
    beta = object$beta

    intercept = object$intercept
    
    Y_hat = X_arr %*% beta + intercept
    if(subtract_mean) {
      Y_hat = Y_hat - mean(Y_hat)
    }
    
    plot(X.seq, Y_hat, type = "l", xlab = "X", ylab = "Y", main = paste0("Spline Fitting for Group", group_index), ...)
}