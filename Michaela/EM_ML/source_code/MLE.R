## Complete data MLE estimators 

mle.mu <- function(x, w){
  sum(x*w)/(sum(w))
}

mle.sigma2 <- function(x, w, mu_mle, nu){
  sum(w*(x-mu_mle)^2) / (nu * length(x))
}
