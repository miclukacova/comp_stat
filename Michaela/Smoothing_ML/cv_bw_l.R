cv_bw_l <- function(x, k = 10, h = seq(0.3, 5, by = 0.1)) {
  
  # Partition
  n <- length(x)
  groups <- sample(rep_len(seq(1,k), length.out = n), replace = FALSE)
  
  # number of obs in each group
  N <- sapply(seq(1,k), function(x) sum(groups == x), simplify = TRUE)
  
  # bandwidth log likelihoods 
  bw_l <- numeric(length(h))
  
  # the kernel density estimates
  f_hat_i <- numeric(n)
  
  for(j in seq_along(h)){
    # Calculating CV density estimates
    for (i in seq_along(x)) {
      k <- groups[i]
      x_m_i <- x[groups != k]
      condition <- (x_m_i - h[j] <= x[i]) & (x[i] <= x_m_i + h[j])
      f_hat_i[i] <- 1 / N[k] * sum((1 - (x[i] - x_m_i)^2 / h[j]^2) * condition)
    }
    f_hat_i <- 3 * f_hat_i / (4 * h[j])
    
    # Calculating log likelihood
    bw_l[j] <- sum(log(f_hat_i))
  }
  
  # Returning the bandwidth that has the largest likelihood 
  return(h[which.max(bw_l)])
}