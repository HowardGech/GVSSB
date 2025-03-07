#' Gives the coverage of a fitted GVSSB model.
#' 
#' This function takes in a fitted GVSSB model and a confidence level, and gives the credible region volume of the fitted model. If a coefficient beta is provided, the coverage of the fitted model is also given. Both group-level and marginal-level credible regions are supported.
#' 
#' @import cluster
#' @param object A fitted GVSSB object.
#' @param alpha A numeric value between 0 and 1 indicating the confidence level. Default is 0.95.
#' @param beta A numeric vector of coefficients. Default is NULL.
#' @param type A character string indicating the type of credible region volume. Default is "group".
#' @return A list of the region volume and coverage (if beta is provided).
#' 
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
#'
#' Y <- X %*% beta + rnorm(n, 0, sd = 3)
#'
#' # fit GVSSB model
#' model.Laplace <- GVSSB(X, Y, groups, prior = 'Laplace')
#' 
#' # get the coverage
#' coverage.group = coverage(model.Laplace, alpha = 0.95, beta = beta, type = "group")
#' coverage.marginal = coverage(model.Laplace, alpha = 0.95, beta = beta, type = "marginal")
#' @export
coverage = function(object, alpha = 0.95, beta = NULL, type = c("group", "marginal")){
  
  # consistency check
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha should be a numeric value between 0 and 1!")
  type = match.arg(type)
  
  # get the fitted model
  groups = object$data$groups
  uni.group = unique(groups)
  G = length(uni.group)
  p_i = sapply(1:G,function(i) sum(groups == uni.group[i]))
  p = length(groups)
  groups_index = sapply(1:p,function(j) which(groups[j]==uni.group))
  gamma = object$gamma
  mu = object$beta
  Sigma = object$Sigma

  if(type == "group"){
    # get the group-level credible region
    
    regionvol = rep(NA, G)
    coverindex = rep(NA, G)
    names(regionvol) = uni.group
    names(coverindex) = uni.group
    for(i in 1:G){
        index = which(groups == uni.group[i])
        # get the region volume
        regionvol[i] = group_volume(gamma[i], 1-alpha, mu[index], Sigma[index,index])
        # get the coverage
        if(!is.null(beta)){
            coverindex[i] = group_coverage(gamma[i], 1-alpha, mu[index], Sigma[index,index], beta[index])
        }
    }
    }

    if(type == "marginal"){
    # get the marginal-level credible region
    regionvol = rep(NA, p)
    coverindex = rep(NA, p)
    for(i in 1:p){
        gamma_index = groups_index[i]
        # get the interval length
        regionvol[i] = marginal_volume(gamma[gamma_index], 1-alpha, mu[i], Sigma[i,i])
        # get the coverage
        if(!is.null(beta)){
            coverindex[i] = marginal_coverage(gamma[gamma_index], 1-alpha, mu[i], Sigma[i,i], beta[i])
        }

    }
    }
    if(!is.null(beta)){
        return(list(regionvol = regionvol, coverindex = coverindex))
    }else{
        return(list(regionvol = regionvol))
    }
}
