GVSSB = function(X, Y, groups, w, sigma_tilde,prior='Gaussian', lambda=1, em=T,alpha= 0,beta = 0,nu=1,threshold = 0.5,tol = 1e-5, iter.max = 5000){
  Y = as.matrix(Y)
  n = nrow(Y)
  Yold = Y
  Y = Y - mean(Y)
  if(!(prior %in% c('Gaussian',"Laplace",'T'))) stop("Improper prior!")
  
  uni.group = unique(groups)

  G = length(uni.group)
  p_i = sapply(1:G,function(i) sum(groups == uni.group[i]))
  p = length(groups)
  groups_index = sapply(1:p,function(x) which(groups[x]==uni.group))
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
  
  if(missing(w)) w = 1/G
  
  ##Use ridge regression as an initializer.
  cvfit = cv.glmnet(X, Y, type.measure = "mse", nfolds = 10, alpha = 0)
  gamma = rep(w, G)
  
  ##discard the intercept.
  fit = glmnet(X, Y, alpha = 0, lambda = cvfit$lambda.min)
  mu = as.numeric(coef(cvfit, cvfit$lambda.min))[-1]
  
  if(missing(sigma_tilde)) sigma_tilde = var(Y[,1]) / qchisq(0.9, df = 4)

  
  XTX = matrix(0, nrow = p, ncol = p)
  for(i in 1:G){
    index = which(groups == uni.group[i])
    XTX[index, index] = t(X[, index]) %*% X[,index]
  }

  
  Sigma = diag(1, p)
  if(prior != 'Gaussian'){
    for(i in 1:G){
      index = groups == uni.group[i]
      Sigma[index, index] = solve(XTX[index, index]/sigma_tilde^2+diag(lambda^2/p_i[i], p_i[i]))
    }
  }
  count = 0
  sigma_old = sigma_tilde
  r = Y - X %*% (gamma[groups_index] * mu)
  while(TRUE){
    count = count+1
    diag_term = rep(1 / lambda^2,G)
    if(prior == 'Laplace') {
      kappa = lambda/sqrt(sapply(1:G, function(n) sum(diag(Sigma[groups == uni.group[n],groups == uni.group[n]]))+sum(mu[groups == uni.group[n]]^2)))
      diag_term = kappa
    }
    if(prior == 'T'){
      kappa = nu*lambda^2+sapply(1:G, function(n) sum(diag(Sigma[groups == uni.group[n], groups == uni.group[n]]))+sum(mu[groups == uni.group[n]]^2))
      diag_term = (p_i + nu) / kappa
    }
    gamma_old = gamma
    mu_old = mu

    
    a = order(sapply(1:G, function(m) sum(mu[groups == uni.group[m]]^2)), decreasing = T)

    
    for(i in a){
      index = which(groups == uni.group[i])
      r = r + X[,index] %*% mu[index] * gamma[i]

      term0 = XTX[index, index] / sigma_tilde^2 + diag(diag_term[i], p_i[i])

      Sigma_m = solve(term0)
      Sigma[index, index] = Sigma_m
      X_m = X[, index]

      mu_m = Sigma_m %*% t(X[,index]) %*% r / sigma_tilde^2
      mu[index] = mu_m
      if(prior == 'Gaussian') logit_gamma = as.numeric(log(w / (1 - w)) + 1/2 * log(det(Sigma_m)) - p_i[i] * log(lambda) + 1/2 * t(mu_m) %*% term0 %*% mu_m)
      if(prior == 'Laplace') logit_gamma = as.numeric((log(det(Sigma_m)) + log(pi)+ p_i[i] * log(lambda^2 / 2) - 2 * log(gamma((p_i[i] + 1) / 2)) - lambda^2  /kappa[i] + t(mu_m) %*% term0 %*% mu_m) / 2 + log(w/ (1 - w)))
      if(prior == 'T') logit_gamma = as.numeric((t(mu_m) %*% term0 %*% mu_m + diag_term[i] * (sum(diag(Sigma_m)) + sum(mu_m^2)) + log(det(Sigma_m)) + 2 * log(gamma((p_i[i] + nu) / 2)) - 2 * log(gamma(nu / 2)) - (p_i[i] + nu) * log(kappa[i] / 2) + nu * log(nu * lambda^2 / 2)) / 2 + log(w/ (1 - w)))
      if(logit_gamma > 100){
        gamma_m = 1
      }else{
        gamma_m = logit_inv(logit_gamma)
      }

      gamma[i] = gamma_m
      r = r - X[,index] %*% mu[index] * gamma[i]
    }
    if(count %% 10 == 0){
      if(em){
        w = mean(gamma)
        lambda_old = lambda
        
        if(prior == 'Gaussian') lambda = sqrt(sum(gamma * sapply(1:G,function(i) sum(diag(Sigma[groups == uni.group[i],groups == uni.group[i]]))+sum(mu[groups == uni.group[i]]^2)))/sum(gamma*p_i))
        if(prior == 'Laplace') lambda = sqrt(sum(gamma*(p_i+1))/sum(gamma*(1/lambda^2+1/kappa)))
        if(prior == 'T') lambda = sqrt(sum(gamma)/sum(gamma*diag_term))
      }
    }

    delta1 = abs(sapply(gamma_old, H) - sapply(gamma, H))
    delta3 = sqrt(mean((mu_old - mu)^2))
    
    if(count %% 10 == 0 & count <= 50){
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
      v = term1+term2+term3
      sigma_tilde = as.numeric(sqrt((v/2 + beta)/(n/2 + alpha)))
      if(abs(sigma_tilde-sigma_old) < 1e-4){
        message(paste("The total number of iterations is:",count,sep=' '))
        break
      }
    }
    if(count %%100 == 0)
      message(paste("The number of current iterations is:",count,sep=' '))
      if(count == iter.max){
        warning("The algorithm doesn't converge!")
        break
      }

  }
  
  
  beta = (gamma[groups_index] > threshold) * mu
  
  betaSD = rep(NA, length(mu))
  for (g in 1 : G) {
    active = which(groups == uni.group[g])
    if (length(active) == 1) {
      betaSD[active] = beta[active] * (sqrt(dim(Xstar)[1]) / sqrt(sum(Xstar[,active]^2)))
    } else {
      betaSD[active] = (Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]])) %*% beta[active]) 
    }
  }
  interceptSD = mean(Yold-Xold%*%betaSD)
  
  return(list(b0=interceptSD, Sigma = Sigma, sigma = sigma_tilde,gamma = gamma, beta=betaSD,lambda=lambda))
}
