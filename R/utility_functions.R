#' @import cluster
#' @importFrom stats coef qchisq qnorm predict
H = function(gamma){
  if(abs(gamma - 1) < 1e-5 | gamma < 1e-5){
    result = 0
  }else{
    result = -gamma * log(gamma) - (1 - gamma) * log(1 - gamma)
  }
  result
}

logit_inv = function(x){
  exp(x) / (1 + exp(x))
}

loss_func = function(loss, vec1, vec2){
  if(!(loss %in% c('L2','L1'))) stop('Loss type incorrect. Only L1 loss and L2 loss are supported.')
  if(loss == 'L1') return(sum(abs(vec1 - vec2)))
  if(loss == 'L2') return(sqrt(sum((vec1 - vec2)^2)))
}

group_coverage = function(gamma,alpha,mu,Sigma,beta){
    if(gamma<=alpha) return(prod(beta==0))
    a.g <- 1 - (1-alpha)/gamma
    if(gamma>=1-alpha){
	    if((mu)%*%solve(Sigma)%*%(mu)<=qchisq(1-a.g, df = length(mu))) {
        if(prod(beta==0)) return(1)
        threshold = qchisq(1 - a.g - (1-gamma), df = length(mu))
        return((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=threshold)
      }
      return((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=qchisq(1-a.g, df = length(mu)))
    }
    if(prod(beta==0)) return(1)
    if((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=qchisq(1 - alpha - (1-gamma), df = length(mu))) return(1)
    return(0)
}

group_volume = function(gamma,alpha,mu,Sigma){
    if(gamma<=alpha) return(0)
    a.g <- 1 - (1-alpha)/gamma
    if(gamma>=1-alpha){
      
	    if((mu)%*%solve(Sigma)%*%(mu)<=qchisq(1-a.g, df = length(mu))) {
        threshold = qchisq(1 - a.g - (1-gamma), df = length(mu))
        ellipsoid = structure(list(cov = Sigma, loc = mu, d2 = threshold),
                   class = "ellipsoid")
    return(volume(ellipsoid))
      }
      ellipsoid = structure(list(cov = Sigma, loc = mu, d2 = qchisq(1-a.g, df = length(mu))),
                   class = "ellipsoid")
    return(volume(ellipsoid))

    }
    ellipsoid = structure(list(cov = Sigma, loc = mu, d2 = qchisq(1 - alpha - (1-gamma), df = length(mu))),
                   class = "ellipsoid")
    return(volume(ellipsoid))
}

marginal_coverage = function(gamma,alpha,mu,Sigma,beta){
    if(gamma<=alpha) return(prod(beta==0))
    a.g <- 1 - (1-alpha)/gamma
    if(gamma>=1-alpha){
	    if((mu)%*%solve(Sigma)%*%(mu)<=qchisq(1-a.g, df = 1)) {
        if(prod(beta==0)) return(1)
        threshold = qchisq(1 - a.g - (1-gamma), df = 1)
        return((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=threshold)
      }
      return((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=qchisq(1-a.g, df = 1))

    }
    if(prod(beta==0)) return(1)
    if((beta-mu)%*%solve(Sigma)%*%(beta-mu)<=qchisq(1 - alpha - (1-gamma), df = 1)) return(1)
    return(0)
}

marginal_volume = function(gamma,alpha,mu,Sigma){
    if(gamma<=alpha) return(0)
    a.g <- 1 - (1-alpha)/gamma
    if(gamma>=1-alpha){
	    if((mu)%*%solve(Sigma)%*%(mu)<=qchisq(1-a.g, df = 1)) {
    return(2*sqrt(Sigma)*qnorm(1 - a.g/2 - (1-gamma)/2))
      }
      return(2*sqrt(Sigma)*qnorm(1-a.g/2))

    }
    return(2*sqrt(Sigma)*qnorm(1 - alpha/2 - (1-gamma)/2))
}