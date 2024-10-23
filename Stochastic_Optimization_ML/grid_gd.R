##################### SGD grid #####################

gd_grid <- function(
    par,
    t0 = 1e-2,
    maxit = 1200,
    cb = NULL,
    epsilon = 1e-5,
    beta = 0.8,
    alpha = 0.1,
    x,
    y,
    ...) {
  
  x_vals <- unique(x)
  matches <- match(x, x_vals)
  n <- length(x)

  for (i in 1:maxit) {
    
    # Computing 
    fs <- f(par, x_vals)[matches]
    nabla_fs <- sapply(seq_along(x_vals), function(i) nabla_f(par, x_vals[i]))
    
    # Calculations of objective and gradient
    value <- sum((y - fs)^2) 
    gr <- - 2 / n * nabla_fs[,matches] %*% (y - fs)
    
    grad_norm <- sum(gr^2)
    
    # Callback
    if (!is.null(cb)) cb()
    
    t <- t0
    # Proposed descent step
    par_new <- par - t * gr
    new_fs <- f(par_new, x_vals)[matches]
    
    # Convergence criterion based on gradient norm
    if (all(abs(par_new - par) <= epsilon)) break
    
    # Backtracking line search
    while (sum((y - new_fs)^2) > value - alpha * t * grad_norm) {
      t <- beta * t
      par_new <- par - t * gr
      new_fs <- f(par_new, x_vals)[matches]
    }
    par <- par_new
  }
  
  if (i == maxit)  warning("Maximal number, ", maxit, ", of iterations reached")
  
  par
}


sgd_grid <- function(
    par,
    grad, # Function of parameter and observation index
    gamma, # Decay schedule or a fixed learning rate
    maxiter = 150, # Max epoch iterations
    sampler = sample, # How data is resampled. Default is a random permutation
    cb = NULL,
    epoch = vanilla,
    m = 1, # Batch size
    x,
    y,
    ...) {
  
  n <- length(x)
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter)
  
  for (k in 1:maxiter) {
    
    if (!is.null(cb)) cb()
    samp <- sampler(n)
    for (j in 1:n) {
      i <- samp[j]
      par <- par - gamma * grad(par, x[i], y[i])
    }
    
  }
  par
}


