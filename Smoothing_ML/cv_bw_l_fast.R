#-------------------------------------------------------------------------------
# CV based optimal band width finder
#-------------------------------------------------------------------------------

cv_bw_l_fast <- function(x, k = 5, h = seq(0.1, 2, by = 0.5)) {
  
  # Partition
  n <- length(x)
  groups <- sample(rep_len(seq(1,k), length.out = n), replace = FALSE)
  
  # number of obs in each group
  N <- sapply(seq(1,k), function(x) sum(groups == x), simplify = TRUE)
  
  # bandwidth log likelihoods and kernel density estimates
  bw_l <- numeric(length(h))
  f_hat_i <- numeric(n)
  
  # term that does not depend on h
  kern_vals <- outer(x, x, function(x1,y1) (x1 - y1)^2)
  
  for(j in seq_along(h)){
    kern_vals_h <- 1 - kern_vals / h[j]^2
    
    # Calculating CV density estimates
    for (i in seq_along(x)) {
      kk <- groups[i]
      condition <- (abs(x - x[i]) <=  h[j]) & (groups != kk)
      f_hat_i[i] <- 1 / N[kk] * sum(kern_vals_h[i, condition])
    }
    f_hat_i <- 3 * f_hat_i / (4 * h[j])
    
    if(any(f_hat_i < 0)){
      warning("Loglikelihood is negative")
    }
    
    # Calculating log likelihood
    bw_l[j] <- sum(log(f_hat_i))
  }
  
  # Returning the bandwidth that has the largest likelihood 
  return(h[which.max(bw_l)])
}