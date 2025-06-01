#' Fit a group linear regression model with discrete spike and slab prior in variational Bayes.
#'
#' This function gives a sparse estimation of the coefficients in group linear model using spike and slab variational Bayes. To make prediction, please use \code{\link{predict.GVSSB}}.
#'
#' @import glmnet
#' @import mvtnorm
#' @import grpreg
#' @importFrom grpreg cv.grpreg
#' @param X An n by p matrix of covariates.
#' @param Y The outcome vector with length n, or the n by 1 outcome matrix.
#' @param groups The group indicator vector with length p.
#' @param w The hyperparameter controlling the sparsity. If missing, the algorithm will use group Lasso to initialize this value.
#' @param sigma_noise The initial value of noise level \eqn{\tilde{\sigma}^2}.
#' @param Omega The prior precision matrix of the coefficients. Default is identity matrix.
#' @param prior The slab part of the coefficient prior. Default is "Gaussian".
#' @param nu Degree of freedom when the slab is T distribution. Only useful when prior='T'. Default is 1.
#' @param lambda The initial value of hyperparameter \eqn{\lambda}. If missing, the algorithm will use group Lasso to initialize this value.
#' @param em Whether update hyperparameters using variational EM algorithm. Default is TRUE.
#' @param alpha Shape of the inverse gamma distribution for the noise level \eqn{\tilde{\sigma}^2}. Default is 0.
#' @param beta Scale of the inverse gamma distribution for the noise level \eqn{\tilde{\sigma}^2}. Default is 0.
#' @param threshold Threshold of the posterior inclusion probability. Only coefficients with posterior inclusion larger than this threshold are considered nonzero (signals).
#' @param tol Convergence tolerance. The algorithm iterates untill the difference of binary entropy between two consecutive steps is smaller than this value. Default is 1e-5.
#' @param iter.min Minimum iteration number of the algorithm. Default is 50.
#' @param iter.max Maximum iteration number of the algorithm. Default is 5000.
#' @param emGap The gap of EM update. Default is 10.
#' @param sigmaGap The gap of noise level update. Default is 10.
#' @param info Whether print the iteration numbers per 100 steps and when the convergence is satisfied. Default is TRUE.
#' @param standardize Whether standardize the covariates. Default is TRUE.
#' @return An object with S3 class "GVSSB".
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
#' @seealso \code{\link{AMSSB}}, \code{\link{predict.GVSSB}}, \code{\link{cv.GVSSB}}
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
#' model.Gaussian <- GVSSB(X, Y, groups, prior = 'Gaussian')
#' model.Laplace <- GVSSB(X, Y, groups, prior = 'Laplace')
#' model.Cauchy <- GVSSB(X, Y, groups, prior = 'T', nu = 1)
#' @export
GVSSB = function(X, Y, groups, Omega, w, lambda, sigma_noise, prior=c('Gaussian','Laplace','T'), nu=1, em=T,alpha= 0,beta = 0, threshold = 0.5,tol = 1e-5, iter.min=50, iter.max = 5000, emGap=10, sigmaGap=10, info = T, standardize=T){
  Y = as.matrix(Y)
  n = nrow(Y)
  if(nrow(X) != n) stop("input dimensions don't match!")
  prior = match.arg(prior)
  uni.group = unique(groups)

  G = length(uni.group)
  p_i = sapply(1:G,function(i) sum(groups == uni.group[i]))
  p = length(groups)
  if(missing(Omega)) Omega = diag(p)
  groups_index = sapply(1:p,function(j) which(groups[j]==uni.group))
  ## Standardize the raw data.
  Yold = Y
  Y = Y - mean(Y)
    if(standardize){
    Xstar = matrix(NA, dim(X)[1], dim(X)[2])
    for (j in 1 : dim(X)[2]) {
        Xstar[,j] = (X[,j] - mean(X[,j]))
    }
    Xtilde = matrix(NA, dim(X)[1], dim(X)[2])
    Qmat = list()
    Dvec = list()

    for (g in 1 : G) {
        active = which(groups == uni.group[g])

        if (length(active) == 1) {
        Xtilde[,active] = sqrt(dim(Xstar)[1]) * (Xstar[,active] / sqrt(sum(Xstar[,active]^2)))
        Qmat[[g]] = NULL
        Dvec[[g]] = NULL
        } else {

        tempX = Xstar[,active]
        SVD = svd((1/n) * t(tempX) %*% tempX)
        Qmat[[g]] = SVD$u
        Dvec[[g]] = SVD$d

        Xtilde[,active] = Xstar[,active] %*% Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]]))
        }
    }
    Xold = X
    X = Xtilde
    }
  Xold = X


  # Use ridge regression as an initializer.
  cvfit_ridge = cv.glmnet(X, Y, type.measure = "mse", nfolds = 10, alpha = 0)
  mu = as.numeric(coef(cvfit_ridge, cvfit_ridge$lambda.min))[-1]

  ## Initialize the hyperparameters if not provided.
  if(missing(lambda) || missing(w)){
    y = Y
    if(info) message('Initializing hyperparameters with group Lasso ...')
      cv_grp = cv.grpreg(X,y,groups,nlambda=10,nfolds=5)
which_lambda_min = which.min(cv_grp$cve)
cv_e_se = cv_grp$cve[which_lambda_min] + cv_grp$cvse[which_lambda_min]
lambda_1se = cv_grp$lambda[cv_grp$cve<=cv_e_se] |> max()
best_model = grpreg(X,y,groups,lambda = lambda_1se)
mu_lasso = as.vector(best_model$beta[-1])
nonzero_index = which(sapply(1:G,function(i) sum(mu_lasso[groups == uni.group[i]]^2)>0))
      if(missing(lambda)){
        if(sum(mu_lasso != 0) != 0){
          if(prior == 'Gaussian') lambda = sqrt((sum(mu_lasso[mu_lasso!=0]%*%Omega[mu_lasso!=0,mu_lasso!=0]%*%mu_lasso[mu_lasso!=0]))/sum(mu_lasso!=0))
          if(prior == 'Laplace') lambda = (sum(mu_lasso!=0)/sum(sapply(1:G,function(i) sqrt(mu_lasso[groups == uni.group[i]]%*%Omega[groups == uni.group[i],groups == uni.group[i]]%*%mu_lasso[groups == uni.group[i]]))))
          if(prior == 'T') lambda=sqrt(sum(sapply(nonzero_index, function(i) p_i[i]+nu*(mu_lasso[groups == uni.group[i]]%*%Omega[groups == uni.group[i],groups == uni.group[i]]%*%mu_lasso[groups == uni.group[i]])))/nu/sum(mu_lasso!=0))
        }
        else{
          if(prior == 'Gaussian') lambda = sqrt((sum(mu%*%Omega%*%mu))/length(mu))
          if(prior == 'Laplace') lambda = (p/sum(sapply(1:G,function(i) sqrt(mu[groups == uni.group[i]]%*%Omega[groups == uni.group[i],groups == uni.group[i]]%*%mu[groups == uni.group[i]]))))
          if(prior == 'T') lambda=sqrt(sum(sapply(1:G, function(i) p_i[i]+nu*(mu[groups == uni.group[i]]%*%Omega[groups == uni.group[i],groups == uni.group[i]]%*%mu[groups == uni.group[i]])))/nu/sum(mu!=0))
        }
      }
      if(missing(w))w = min(max(mean(sapply(1:G, function(i) sum(mu_lasso[groups==uni.group[i]]^2)>0)),1/G),1-1/G)
  }
  gamma = rep(w, G)
  if(missing(sigma_noise)) sigma_noise = mean((Y - X %*% mu)^2)
  sigma_tilde = sigma_noise
  ## Avoid duplicate computation.
  XTX = matrix(0, nrow = p, ncol = p)
  for(i in 1:G){
    index = which(groups == uni.group[i])
    XTX[index, index] = t(X[, index]) %*% X[,index]
  }

  ## Initialize variational covariance matrix.
  Sigma = diag(1, p)

    for(i in 1:G){
      index = groups == uni.group[i]
      if(prior == 'Gaussian') Sigma[index, index] = solve(XTX[index, index]/sigma_tilde^2+1/lambda^2*Omega[index,index])
      if(prior == 'Laplace') Sigma[index, index] = solve(XTX[index, index]/sigma_tilde^2+lambda^2*2*Omega[index,index])
      if(prior == 'T') Sigma[index, index] = solve(XTX[index, index]/sigma_tilde^2+1/lambda^2*Omega[index,index])
    }

  count = 0
  sigma_old = sigma_tilde

  kappa = rep(1,G)
  r = Y - X %*% (gamma[groups_index] * mu)
  while(TRUE){
    count = count+1

    ## Update kappa for Laplace and T distribution.
    if(prior == 'Gaussian') diag_term = rep(1 / lambda^2,G)
    else {
      kappa = sapply(1:G, function(n) sum(diag(Omega[groups == uni.group[n],groups == uni.group[n]] %*% Sigma[groups == uni.group[n],groups == uni.group[n]]))+sum(mu[groups == uni.group[n]]%*%Omega[groups == uni.group[n],groups == uni.group[n]]%*%mu[groups == uni.group[n]]))
    }
    if(prior == 'Laplace') {
      diag_term = lambda/sqrt(kappa)
    }
    if(prior == 'T'){
      diag_term = (p_i + nu) / (kappa + lambda^2 * nu)
    }
    gamma_old = gamma
    mu_old = mu

    ## Determine the optimization order for coordinate ascent.
    a = order(sapply(1:G, function(m) sum(mu[groups == uni.group[m]]^2)), decreasing = T)

    ## Main part.
    for(i in a){
      index = which(groups == uni.group[i])
      X_m = matrix(X[, index],n)

      r = r + X_m %*% mu[index] * gamma[i]

      Sigma_m_inv = XTX[index, index] / sigma_tilde^2 + diag_term[i]*Omega[index, index]

      Sigma_m = solve(Sigma_m_inv)
      Sigma[index, index] = Sigma_m


      mu_m = Sigma_m %*% t(X_m) %*% r / sigma_tilde^2
      mu[index] = mu_m
      if(prior == 'Gaussian') logit_gamma = as.numeric(log(w / (1 - w)) + 1/2 * log(det(Omega[index,index]%*%Sigma_m)) - p_i[i] * log(lambda) + 1/2 * t(mu_m) %*% Sigma_m_inv %*% mu_m)
      if(prior == 'Laplace') logit_gamma = as.numeric((log(det(Omega[index,index]%*%Sigma_m)) + log(3.14159)+ p_i[i] * log(lambda^2 / 2) - 2 * log(gamma((p_i[i] + 1) / 2)) - lambda*sqrt(kappa[i]) + t(mu_m) %*% Sigma_m_inv %*% mu_m) / 2 + log(w/ (1 - w)))
      if(prior == 'T') logit_gamma = as.numeric((t(mu_m) %*% Sigma_m_inv %*% mu_m + log(det(Omega[index,index]%*%Sigma_m)) + 2 * log(gamma((p_i[i] + nu) / 2)) - 2 * log(gamma(nu / 2)) - (p_i[i] + nu) * log((nu*lambda^2+kappa[i]) / 2) + nu * log(nu * lambda^2 / 2) + kappa[i]*diag_term[i]) / 2 + log(w/ (1 - w)))
      if(logit_gamma > 100){
        gamma_m = 1
      }else{
        gamma_m = logit_inv(logit_gamma)
      }
      gamma[i] = gamma_m
      r = r - X_m %*% mu[index] * gamma[i]
    }

    ## Variational EM update for hyperparameters.
    if(count %% emGap == 0){
      if(em){
        w = mean(gamma)
        lambda_old = lambda

        if(prior == 'Gaussian') lambda = sqrt(sum(gamma * sapply(1:G,function(i) sum(diag(Omega[groups == uni.group[i],groups == uni.group[i]]%*%Sigma[groups == uni.group[i],groups == uni.group[i]])) + sum(mu[groups == uni.group[i]]%*%Omega[groups == uni.group[i],groups == uni.group[i]]%*%mu[groups == uni.group[i]]))) / sum(gamma * p_i))
        if(prior == 'Laplace') lambda = sqrt(sum(gamma * (p_i + 1)) / sum(gamma*(sqrt(kappa)/lambda+1/lambda^2)))
        if(prior == 'T') lambda = sqrt(sum(gamma) / sum(gamma * diag_term))
      }
    }

    delta1 = abs(sapply(gamma_old, H) - sapply(gamma, H))
    delta3 = sqrt(mean((mu_old - mu)^2))


    ## Update noise level.
    if(count %% sigmaGap == 0 & count <= iter.min){
      sigma_old = sigma_tilde
      term1 = sum(r^2)
      term2 = sum(sapply(1:G, function(n) gamma[n] * (1 - gamma[n]) * t(mu[groups == uni.group[n]]) %*% XTX[groups == uni.group[n],groups == uni.group[n]] %*% mu[groups == uni.group[n]]))
      term3 = sum(sapply(1:G, function(n) gamma[n] * sum(diag(XTX[groups == uni.group[n],groups == uni.group[n]] %*% Sigma[groups == uni.group[n],groups == uni.group[n]]))))
      v = term1+term2+term3
      sigma_tilde = as.numeric(sqrt((v/2 + beta)/(n/2 + alpha)))
    }
    if(count > 50 & max(c(delta1)) < tol){
      sigma_old = sigma_tilde
      term1 = sum(r^2)
      term2 = sum(sapply(1:G, function(n) gamma[n] * (1 - gamma[n]) * t(mu[groups == uni.group[n]]) %*% XTX[groups == uni.group[n],groups == uni.group[n]] %*% mu[groups == uni.group[n]]))
      term3 = sum(sapply(1:G, function(n) gamma[n] * sum(diag(XTX[groups == uni.group[n],groups == uni.group[n]] %*% Sigma[groups == uni.group[n],groups == uni.group[n]]))))
      v = term1 + term2 + term3
      sigma_tilde = as.numeric(sqrt((v/2 + beta)/(n/2 + alpha)))
      if(abs(sigma_tilde - sigma_old) < 1e-4){
        if(info) message(paste("The total number of iterations is:",count,sep=' '))
        break
      }
    }
    if(count %% 100 == 0)
      if(info) message(paste("The number of current iterations is:",count,sep=' '))
      if(count == iter.max){
        warning("The algorithm doesn't converge!")
        break
      }

  }

  ## Compute the estimation of coefficients and intercept for raw data.
  beta_raw = rep(NA, length(mu))
  mu_raw = rep(NA, length(mu))
  Sigma_raw = matrix(0, nrow = length(mu), ncol = length(mu))
  if(standardize){
  for (g in 1 : G) {
    active = which(groups == uni.group[g])
    if (length(active) == 1) {
      beta_raw[active] = mu[active] * (sqrt(dim(Xstar)[1]) / sqrt(sum(Xstar[,active]^2)))
      mu_raw[active] = beta_raw[active]
      Sigma_raw[active,active] = Sigma[active, active] * (sqrt(dim(Xstar)[1]) / sqrt(sum(Xstar[,active]^2)))^2
    } else {
      beta_raw[active] = (Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]])) %*% mu[active])
      mu_raw[active] = beta_raw[active]
      Sigma_raw[active,active] = Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]])) %*% Sigma[active, active] %*% diag(1/sqrt(Dvec[[g]])) %*% t(Qmat[[g]])
    }
  }
  }else{
    beta_raw = mu
    mu_raw = mu
    Sigma_raw = Sigma
  }
  beta_raw = beta_raw * (gamma[groups_index] > threshold)
  intercept_raw = mean(Yold-Xold%*%beta_raw)
  selected_groups = which(gamma > threshold)
  data_model = list(Omega = Omega, groups = groups)

  processed_model  = list(data = data_model, intercept=intercept_raw, sigma_noise = sigma_tilde, mu = mu_raw, Sigma = Sigma_raw,  gamma = gamma, beta=beta_raw, selected_groups = selected_groups, lambda=lambda, prior = prior, nu = nu)
  class(processed_model) = 'GVSSB'

  return(processed_model)
}