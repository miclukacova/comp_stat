vec_f_random <- function(m, mu, sigma, f, alpha_p) {
  y <- rnorm(m, mu, sigma)
  u <- runif(m)
  accept <- u <= f(y) / dnorm(y, mean = mu, sd = sigma) * alpha_p
  
  return("Sample" = y[accept])
}

rng_wrapper <- function(rng, fact = 1.05, M_min = 100) {
  function(N, ...) {
    j <- 0     # The number of iterations
    l <- 0     # The number of accepted samples
    counter <- 0 # The number of proposals
    
    x <- list()
    
    while (l < N) {
      j <- j + 1 
      M <- floor(max(fact * (N - l), M_min))
      x[[j]] <- rng(M, ...)
      counter <- counter + M
      l <- l + length(x[[j]])
      # Update 'fact' by estimated acceptance probability l / n
      if (j == 1) fact <- fact * N / l
    }
    return(list("sample" = unlist(x)[1:N], "accept" = length(unlist(x)) / counter))
  }
}
