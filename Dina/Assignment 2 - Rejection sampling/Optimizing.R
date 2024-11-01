## Function factory for rejection sampling
new_rejection_sampler <- function(generator) {
  function(n, ...) {
    alpha <- 1
    y <- numeric(0)
    n_accepted <- 0
    while (n_accepted < n) {
      m <- ceiling((n - n_accepted) / alpha)
      y_new <- generator(m, ...)
      n_accepted <- n_accepted + length(y_new)
      if (length(y) == 0) {
        alpha <- (n_accepted + 1) / (m + 1) # Estimate of alpha
      }
      y <- c(y, y_new)
    }
    list(x = y[seq_len(n)], alpha = alpha)
  }
}


## The parallelized one
# Slopes
a_i <- function(x_i) tar_dens_log_difference(x_i) 

# Intercepts
b_i <- function(x_i, a_i) log(tar_dens(x_i)) - a_i * x_i

# Interval points
z_i <- function(a1, a2, b1, b2) (b2 - b1) / (a1 - a2)

# R_i's
r_i <- function(as, bs, zs, n) {
  1 / as * exp(bs) * (exp(as * zs[2:(n+1)]) - exp(as * zs[1:n]))
}

piece_lin_rejec_samp_par <- function(N, ys) {
  
  # Calculating a's, b's, z's
  as <- sapply(ys, a_i, simplify = TRUE)
  bs <- mapply(FUN = b_i, ys, as)
  n <- length(bs)
  zs <- c(-Inf, mapply(FUN = z_i, as[1:(n-1)], as[2:n], bs[1:(n-1)], bs[2:n]), Inf)
  
  # Bookkeeping
  # I_i integrals
  R <- r_i(as, bs, zs, n)
  
  # Distribution function (ish)
  Q <- c(0, cumsum(R))
  
  # Drawing from piecewise linear density and uniform
  u0 <- Q[n + 1] * runif(N)
  u <- runif(N)
  
  # Determine the interval that each point belongs to
  geq_z <-outer(u0, Q[1:n], FUN = function(y1, y2) y1 > y2)
  leq_z <-outer(u0, Q[2:(n+1)], FUN = function(y1, y2) y1 <= y2)
  I <- geq_z & leq_z
  
  
  x <- numeric(N)
  accept <- logical(N)
  
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Create storage for results
  results <- foreach(i = 1:N, .combine = rbind, .export = c("tar_dens", "I", "u0", "Q", "as", "bs", "zs", "zx", "poisson_data")) %dopar% {
    # Finding the interval x_i belongs to
    int <- which(I[i,] == 1)
    
    # Make sure 'int' is not empty
    if (length(int) == 0) {
      return(c(NA, NA))  # Return NA if no interval found
    }
    
    # Taking the inverse cdf
    x_i <- log((u0[i] - Q[int]) * as[int] * exp(-bs[int]) + exp(as[int] * zs[int])) / as[int]
    
    # Acceptance step
    accept_i <- u[i] <= tar_dens(x_i) / exp(as[int] * x_i + bs[int])
    
    # Return the results as a vector
    c(x_i, accept_i)
  }
  
  # Ensure results have correct dimensions before indexing
  if (ncol(results) == 2) {
    x <- results[, 1]  # First column for x
    accept <- results[, 2]  # Second column for accept
  } else {
    stop("Results do not have the expected structure.")
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Check if x and accept are not NULL or NA
  if (any(is.na(x)) || any(is.na(accept))) {
    warning("Some values in x or accept are NA. Check input data and function implementation.")
  }
  
  
  
  return(as.vector(x[accept]))
}

aff_rej_par <- new_rejection_sampler(piece_lin_rejec_samp_par)


## One with the loop in Rcpp







