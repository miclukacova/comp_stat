
##################################################
################# Gaussian Envelope ##############
##################################################

# Ratio to minimize
ratio <- function(y, mu, sigma) 1 / (sqrt(2 * pi *sigma^2)) * 
  exp(- (y - mu)^2 / ( 2 * sigma^2) - y * sum1  + sum(exp(y * x)))

log_ratio <- function(y, mu, sigma) - log(sqrt(2 * pi *sigma^2)) -
  (y - mu)^2 / ( 2 * sigma^2) - y * sum1  + sum(exp(y * x))
 
f_sample <- function(n) {
  y <- numeric(n)
  counter <- 0
  for (i in 1:n) {
    reject <- TRUE
    while (reject) {
      y0 <- rnorm(1, mu_opt, sigma_opt)
      u <- runif(1)
      reject <- u > f_star(y0) / g_star(y0) * alpha_p
      counter <- counter + 1
    }
    y[i] <- y0
  }
  return(list("Sample" = y, "Accept" = n / counter))
}


# Vectorized
vec_f_random <- function(m) {
  y <- rnorm(m, mu_opt, sigma_opt)
  u <- runif(m)
  y1 <- sapply(y, f_star, simplify = TRUE)
  accept <- u <= y1 / g_star(y) * alpha_p
  
  return("Sample" = y[accept])
}

# Wrapper function
rng_wrapper <- function(rng, fact = 1.05, M_min = 100) {
  function(N, ...) {
    j <- 0     # The number of iterations
    l <- 0     # The number of accepted samples
    counter <- 0 # The number of proposals
    
    x <- list()
    
    while (l < N) {
      j <- j + 1 # Adapt the number of proposals
      M <- floor(max(fact * (N - l), M_min))
      x[[j]] <- rng(M, ...)
      counter <- counter + M
      l <- l + length(x[[j]])
      # Update 'fact' by estimated acceptance probability l / n
      if (j == 1) fact <- fact * N / l
    }
    return(list("Sample" = unlist(x)[1:N], "Accept" = length(unlist(x)) / counter))
  }
}

