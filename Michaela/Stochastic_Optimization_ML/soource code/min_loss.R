# Minimal loss

min_loss <- function(x, y, decay = decay_scheduler, sv){
  decayy <- decay_scheduler(gamma0 = 1, a = 0.6, n1 = 100, gamma1 = 0.001)
  
  optim_sgd <- SGD(par0 = sv, grad = grad_rcpp, 
                   gamma = decayy, x = x, y = y, epsilon = 1e-8, maxit = 400)
  
  loss <- H_mult(x = x, y = y, 
                 alpha = optim_sgd$trace$par.1, 
                 beta = optim_sgd$trace$par.2,
                 gamma = optim_sgd$trace$par.3,
                 rho = optim_sgd$trace$par.4)
  
  min_loss <- loss %>% min()
  optim_par <- optim_sgd$trace[which(min_loss == loss), 1:4]
  
  return(list(min_loss, optim_par))
}



