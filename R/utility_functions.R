H = function(gamma){
  if(abs(gamma-1)<1e-5 | gamma<1e-5){
    result = 0
  }else{
    result = -gamma*log(gamma)-(1-gamma)*log(1-gamma)
  }
  result
}

logit_inv = function(x){
  exp(x)/(1+exp(x))
}
