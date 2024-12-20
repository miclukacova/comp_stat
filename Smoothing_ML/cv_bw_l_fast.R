#-------------------------------------------------------------------------------
# CV based optimal band width finder
#-------------------------------------------------------------------------------

cv_bw_l_fast <- function(x, 
                         k = 10,
                         h = seq(0.3, 5, by = 0.1), 
                         loop = inner_loop,
                         do_parallel = FALSE) {
  
  # Partition
  n <- length(x)
  groups <- sample(rep_len(seq(1,k), length.out = n), replace = FALSE)
  
  # number of obs in each group
  N <- sapply(seq(1,k), function(x) sum(groups == x), simplify = TRUE)
  
  # bandwidth log likelihoods and kernel density estimates
  bw_l <- numeric(length(h))
  
  # term that does not depend on h
  kern_vals <- outer(x, x, function(x1,y1) (x1 - y1)^2)
  
  if(do_parallel == TRUE){
    
    bw_l_results <- mclapply(seq_along(h), function(j) {
      kern_vals_h <- 1 - kern_vals / h[j]^2
      
      # Calculating CV density estimates
      f_hat_i <- loop(x, h[j], groups, kern_vals_h, N)
      
      # Calculating log likelihood for bandwidth h[j]
      log_likelihood <- sum(log(f_hat_i))
      
      # Return log likelihood for this bandwidth
      return(log_likelihood)
    }, mc.cores = detectCores() - 1)
    
    bw_l <- unlist(bw_l_results)
  }
  else{
    for(j in seq_along(h)){
      kern_vals_h <- 1 - kern_vals / h[j]^2
      
      f_hat_i <- inner_loop(x, h[j], groups, kern_vals_h, N)
      
      print(any(f_hat_i == 0))
    
      # Calculating log likelihood
      bw_l[j] <- sum(log(f_hat_i))
    }
  }
  # Returning the bandwidth that has the largest likelihood 
  return(h[which.max(bw_l)])
}


inner_loop <- function(x, h, groups, kern_vals_h, N){
  
  f_hat_i <- numeric(length(x))
  
  # Calculating CV density estimates
  for (i in seq_along(x)) {
    kk <- groups[i]
    condition <- (abs(x - x[i]) <=  h) & (groups != kk)
    f_hat_i[i] <- 1 / N[kk] * sum(kern_vals_h[i, condition])
  }
  f_hat_i <- 3 * f_hat_i / (4 * h)
  
  return(f_hat_i)
}


