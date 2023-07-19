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