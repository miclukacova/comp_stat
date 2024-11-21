
# Vanilla tuning
vanilla_tuning <- function(tuning_rounds = 50, iterations = 250, seed = 12345){
  
  Loss <- Inf
  set.seed(seed)
  
  for (i in 1:tuning_rounds){
    
    gamma1 <- rexp(1, 5)
    decay <- decay_scheduler(gamma0 = 1, a = 0.75, n1 = iterations, gamma1 = gamma1)
    
    vanilla_SGD_object <- SGD(par0 = init_par, 
                            grad = gradient_rcpp, 
                            gamma = decay,
                            epoch = vanilla,
                            maxiter = iterations, 
                            sampler = sample,
                            x = x, y = y,
                            true_par = true_par)
    
    if(Loss > H(x, y, vanilla_SGD_object$est)){
      Loss <- H(x, y, vanilla_SGD_object$est)
      best_vanilla <- vanilla_SGD_object
    }
  }
  return(best_vanilla)
}
#vanilla_tuning(15)



# Batch tuning

batch_tuning <- function(tuning_rounds = 50, iterations = 250, seed = 12345){
  
  Loss <- Inf
  set.seed(seed)
  
  for (i in 1:tuning_rounds){
    
    m <- sample(25:200, 1)
    gamma1 <- rexp(1, 10)
    decay <- decay_scheduler(gamma0 = 1, a = 0.75, n1 = iterations, gamma1 = gamma1)
    
    batch_SGD_object <- SGD(par0 = init_par, 
                            grad = gradient_rcpp, 
                            gamma = decay,
                            epoch = batch, 
                            m = m, 
                            maxiter = iterations, 
                            sampler = sample,
                            x = x, y = y,
                            true_par = true_par)
    
    if(Loss > H(x, y, batch_SGD_object$est)){
      Loss <- H(x, y, batch_SGD_object$est)
      best_batch <- batch_SGD_object
    }
  }
  return(best_batch)
}
#batch_tuning(25)




# Momentum tuning
momentum_tuning <- function(tuning_rounds = 50, iterations = 250, seed = 12345){
  
  Loss <- Inf
  set.seed(seed)
  
  for (i in 1:tuning_rounds){
    
    m <- sample(25:200, 1)
    gamma1 <- rexp(1, 10)
    decay <- decay_scheduler(gamma0 = 1, a = 0.75, n1 = iterations, gamma1 = gamma1)
    beta <- runif(1, 0.8, 0.95)
    
    momentum_SGD_object <- SGD(par0 = init_par, 
                            grad = gradient_rcpp, 
                            gamma = decay,
                            beta = beta,
                            epoch = momentum(), 
                            m = m, 
                            maxiter = iterations, 
                            sampler = sample,
                            x = x, y = y,
                            true_par = true_par)
    
    if(Loss > H(x, y, momentum_SGD_object$est)){
      Loss <- H(x, y, momentum_SGD_object$est)
      best_momentum <- momentum_SGD_object
    }
  }
  return(best_momentum)
}
#momentum_tuning(25)



# Adam tuning
adam_tuning <- function(tuning_rounds = 50, iterations = 250, seed = 12345){
  
  Loss <- Inf
  set.seed(seed)
  
  for (i in 1:tuning_rounds){
    
    m <- sample(100:500, 1)
    gamma1 <- rexp(1, 3)
    decay <- decay_scheduler(gamma0 = 1, a = 0.75, n1 = iterations, gamma1 = gamma1)
    beta1 <- runif(1, 0.8, 0.95)
    beta2 <- runif(1, 0.8, 0.95)
    
    adam_SGD_object <- SGD(par0 = init_par, 
                            grad = gradient_rcpp, 
                            gamma = decay,
                            beta1 = beta1,
                            beta2 = beta2,
                            epoch = adam(), 
                            m = m, 
                            maxiter = iterations, 
                            sampler = sample,
                            x = x, y = y,
                            true_par = true_par)
    
    if(Loss > H(x, y, adam_SGD_object$est)){
      Loss <- H(x, y, adam_SGD_object$est)
      best_adam <- adam_SGD_object
    }
  }
  return(best_adam)
}
#adam_tuning(25)



GD_tuning <- function(tuning_rounds = 50, iterations = 250, seed = 12345){
  
  Loss <- Inf
  set.seed(seed)
  
  for (i in 1:tuning_rounds){
    
    gamma <- rexp(1, 1)
    beta <- runif(1, 0.5, 0.95)
    d <- runif(1, 0.5, 0.95)
    
    GD_object <- GD(par0 = init_par, 
                    H = H_rcpp,
                    grad = gradient_rcpp, 
                    d = d,
                    gamma = gamma,
                    beta = beta,
                    momentum = TRUE,
                    clipping = TRUE,
                    maxiter = iterations,
                    x = x, y = y,
                    true_par = true_par)
    
    if(Loss > H(x, y, GD_object$est)){
      Loss <- H(x, y, GD_object$est)
      best_GD <- GD_object
    }
  }
  return(best_GD)
}
#GD_tuning(25)
