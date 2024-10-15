#------------------------------------------------------------------------------#
#              Different implementations of the gaussian envelope              #
#------------------------------------------------------------------------------#

gaussian_loop <- function(N, mu, sigma, alpha_star){
  # Initialize samples
  samples <- numeric(N)
  
  # Count number of rejections
  rejection_count <- 0
  
  for (i in 1:N){
    reject <- TRUE
    
    while (reject) {
      x0 <- rnorm(1, mean = mu, sd = sigma) # Sample from proposal distribution
      u <- runif(1) # Sample from uniform distribution
      
      target_x0 <- target_distribution_pois(x0)  # Evaluate target distribution at x0
      proposal_x0 <- dnorm(x0, mean = mu, sd = sigma) # Evaluate proposal distribution at x0
      
      reject <- u > alpha_star * target_x0/proposal_x0 # Acceptance-rejection step
      
      if (reject) {
        rejection_count <- rejection_count + 1 # Count rejections
      }
      
    }
    samples[i] <- x0 # Store accepted sample
  }
  alpha_hat <- N / (N + rejection_count) # Estimate acceptance rate
  return(list(samples = samples, rejection_count = rejection_count, alpha_hat = alpha_hat))
}


gaussian_fast_loop <- function(N, mu, sigma, alpha_star){
  # Initialize samples
  samples <- numeric(N)
  
  # Initalize poisson variables
  x_poisson <- poisson$x
  z_poisson <- poisson$z
  xz_poisson <- x_poisson*z_poisson
  
  # Count number of rejections
  rejection_count <- 0
  
  for (i in 1:N){
    reject <- TRUE
    
    while (reject) {
      x0 <- rnorm(1, mean = mu, sd = sigma) # Sample from proposal distribution
      u <- runif(1) # Sample from uniform distribution
      
      target_x0 <- rcpp_target_distribution_pois(x0, x = x_poisson, z = z_poisson, xz = xz_poisson)  # Evaluate target distribution at x0
      proposal_x0 <- dnorm(x0, mean = mu, sd = sigma) # Evaluate proposal distribution at x0
      
      reject <- u > alpha_star * target_x0/proposal_x0 # Acceptance-rejection step
      
      if (reject) {
        rejection_count <- rejection_count + 1 # Count rejections
      }
      
    }
    samples[i] <- x0 # Store accepted sample
  }
  alpha_hat <- N / (N + rejection_count) # Estimate acceptance rate
  return(list(samples = samples, rejection_count = rejection_count, alpha_hat = alpha_hat))
}


gaussian_random <- function(N, mu, sigma, alpha_star){

  # Initalize poisson variables
  x_poisson <- poisson$x
  z_poisson <- poisson$z
  xz_poisson <- x_poisson*z_poisson
  
  x0 <- rnorm(N, mean = mu, sd = sigma) # Sample from proposal distribution
  u <- runif(N) # Sample from uniform distribution
  target_x0 <- rcpp_target_distribution_pois(x0, x = x_poisson, z = z_poisson, xz = xz_poisson)  # Evaluate target distribution at x0
  proposal_x0 <- dnorm(x0, mean = mu, sd = sigma) # Evaluate proposal distribution at x0
  accept <- u <= alpha_star * target_x0/proposal_x0 # Acceptance-rejection step
  accepted <- x0[accept]
  
  #alpha_hat <- length(accepted) / N # Estimate acceptance rate
  return(accepted)
}




rng_vec <- function(rng, fact = 1.2, M_min = 100) {
  CSwR::force_all()
  function(N, ..., cb) {
    j <- 0     # The number of iterations 
    l <- 0     # The number of accepted samples
    x <- list() # The list of accepted samples
    
    total_samples <- 0
    
    while (l < N) {
      j <- j + 1
      
      # Adapt the number of proposals
      M <- floor(max(fact * (N - l), M_min))
      
      x[[j]] <- rng(M, ...)
      
      # Update total number of samples
      total_samples <- total_samples + M
      
      l <- l + length(x[[j]])
      if (!missing(cb)) cb() # Callback
      # Update 'fact' by estimated acceptance probability l / n
      if (j == 1) fact <- fact * N / l
    }
    
    samples <- unlist(x)
    
    return(list(samples = samples[1:N], alpha_hat = length(samples) / total_samples))
  }
}



gaussian_vec <- rng_vec(gaussian_random)


















N <- 1000
break_points <- c(0.1, 0.2, 0.3)
log_target_dist <- log_target_distribution_pois
log_target_dist_diff <- log_target_distribution_diff_pois

adaptive_rejection_sampler <- function(N, break_points, log_target_dist, log_target_dist_diff){
  
  
  # Calculate a's, b's and z's
  a <- log_target_dist_diff(break_points)
  if (sum(duplicated(a)) != 0 || sum(a == 0))
    stop("\n
    The implementation requires the a's to be different and non-zero. 
    Choose different break points to ensure this."
    )
  
  b <- log_target_dist(break_points) - a * break_points
  z <- z_function(a, b, -Inf, Inf)
  
  
  # Use the calculated values of a, b and z to calculate the Q_i's
  R_is <- H_i(a, b, z, break_points)
  Q_is <- cumsum(R_is)
  
  # Return d as the last element of Q_is
  d <- Q_is[length(Q_is)]
  
  # Simulation from the proposal distribution
  
  # Check if dq_j is in a (Q_{i-1}, Q_i] interval and store this information in indicator I with dimension m x N
  dq <- d * runif(N)
  I <- dq_in_Q_i(dq, Q_is, z)
  
  # If dq is in the interval the inverse of G is given by x(q) = H_i^{-1}(dq - Q_{i-1})
  
  #Note dq_matrix contains a lot of 0's (since dq is not in the interval for all i's) and we only need the non-zero elements, so we convert it into a sparse matrix.
  dq_matrix <- as(t(t(I) * dq), "dgCMatrix") # This could be a bottleneck as we transpose the I matrix twice
  dq_matrix_non_zero <- dq_matrix@x # Extract non-zero elements of dq
  dq_matrix_non_zero_rows <- dq_matrix@i + 1 # Extract the row indices of the non-zero elements of dq
  
  # We extract the non-zero elements of dq_matrix and the row indicies of the non-zero elements of dq_matrix, so we know which values of Q_{i-1}, z_{i-1}, a_i and b_i to use for each non-zero element of dq.
  
  # For each non-zero element of dq, we return the corresponding Q_{i-1}, z_{i-1}, a_i and b_i and calculate the inverse of H_i
  Q_is_minus1 <- c(0, Q_is[-length(Q_is)])[dq_matrix_non_zero_rows]
  z_is_minus1 <- z[-length(z)][dq_matrix_non_zero_rows]
  a_is <- a[dq_matrix_non_zero_rows]
  b_is <- b[dq_matrix_non_zero_rows]
  
  # Finally we sample from the proposal distribution
  x <- 1 / a_is * log(a_is * exp(-b_is) * (dq_matrix_non_zero - Q_is_minus1) + exp(a_is*z_is_minus1))
  
  
  accept <- logical(N)
  u <- runif(N)
  
  for (i in 1:length(a)){ #1 to m. Note since m is usually quite small, there is not much to gain by removing this for loop
    group_i <- as.logical(I[i,]) # Index to select observations belonging to interval i
    
    accept[group_i] <- u[group_i] <= exp(log_target_dist(x[group_i]) - a[i] * x[group_i] - b[i])
  }
  
  return(x[accept])
}