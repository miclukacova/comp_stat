###############################################
#                                             #
#       Bandwidth selection algorithms        #
#                                             #
###############################################


AMISE_epa_bandwidth <- function(x){
  n <- length(x)
  K_square_norm <- 3/5
  k_sigma_4 <- 1/25
  f_tilde_norm2 <- epanechnikov_f_tilde_norm2(x)
  
  (K_square_norm/(k_sigma_4 * f_tilde_norm2))^(1/5) * n^(-1/5)
}

AMISE_epa_bandwidth_rcpp <- function(x){
  n <- length(x)
  K_square_norm <- 3/5
  k_sigma_4 <- 1/25
  f_tilde_norm2 <- epanechnikov_f_tilde_norm2_rcpp(x)
  
  (K_square_norm/(k_sigma_4 * f_tilde_norm2))^(1/5) * n^(-1/5)
}








#------------ Cross-validation ------------



#Make cv class?

cv_naive <- function(x, h_grid, folds, kernel = epanechnikov, seed = NA) {
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  # Define length of x
  N <- length(x)
  group <- sample(rep(1:folds, length.out = N))
  
  # List of cross-validated sums to compare values of h in the end
  cv_sum <- numeric(length(h_grid))
  
  
  # Loop over h-values
  for (h_index in seq_along(h_grid)){
    
    # Get value of h from grid of h's
    h <- h_grid[h_index]
    
    # Define vector to store log-likelihood value of each fold
    fold_l_cv <- numeric(folds)
    
    
    for (fold in seq_along(1:folds)){
      
      # Split data into folds based on group partition earlier
      other_sets <- x[group != fold]
      hold_out_set <- x[group == fold]
      
      # Define vector to store f_i_hat for each x_i in hold_out_set
      f_i_hats <- numeric(length(hold_out_set))
      N_i <- length(other_sets)
      
      for (x_i in seq_along(hold_out_set)){
        f_i_hat <- 1/(h * N_i) * sum(epanechnikov((hold_out_set[x_i] - other_sets) / h))
        f_i_hats[x_i] <- f_i_hat
      }
      fold_l_cv[fold] <- sum(log(f_i_hats))
    }
    cv_sum[h_index] <- sum(fold_l_cv)
  }
  
  cv_result <- data.frame(h = h_grid, cv_result = cv_sum)
  best_h <- h_grid[which.max(cv_sum)]
  
  return(list(cv_result = cv_result, best_h = best_h))
}


cv_outer <- function(x, h_grid, folds, kernel = epanechnikov, seed = NA) {
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  # Initiate variables
  group <- sample(rep(1:folds, length.out = length(x))) # Randomly assign each observation to a fold
  ll_sum <- numeric(length(h_grid))                     # Vector to store ll for each h
  
  
  for (h_index in seq_along(h_grid)){
    
    h <- h_grid[h_index]          # Get h from grid of h's
    fold_l_cv <- numeric(folds)   # Vector to store log-likelihood for each fold
    
    for (fold in seq_along(1:folds)) {
      # Define hold-out set and remaining sets for fold
      other_sets <- x[group != fold]
      hold_out_set <- x[group == fold]
      
      N_i <- length(other_sets)
      
      # Precalculate kernel input
      kernel_input <- Rfast::Outer(x = hold_out_set, y = other_sets, "-") / h
      
      # Calculate i'th kernel estimate
      f_i_hats <- colSums(kernel(kernel_input))
      
      
      # Get log-likelihood for fold
      fold_l_cv[fold] <- sum(log(1/(h * N_i) * f_i_hats))
    }
    
    # Sum log-likelihoods to obtain log-likelihood for h
    ll_sum[h_index] <- sum(fold_l_cv)
    
  }
  
  cv_result <- data.frame(h = h_grid, cv_result = ll_sum)
  best_h <- h_grid[which.max(ll_sum)]
  
  return(list(cv_result = cv_result, best_h = best_h))
}



cv_Outer <- function(x, h_grid, folds, kernel = epanechnikov, seed = NA) {
  
  N <- length(x)
  
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  group <- sample(rep(1:folds, length.out = N))
  cv_sum <- numeric(length(h_grid))
  
  outer_input <- Rfast::Outer(x = x, y = x, "-")
  
  for (h_index in seq_along(h_grid)){
    h <- h_grid[h_index]
    fold_l_cv <- numeric(folds)
    
    total_input <- outer_input / h
    
    for (fold in seq_along(1:folds)) {
      other_sets <- x[group != fold]
      hold_out_set <- x[group == fold]
      
      N_i <- length(other_sets)
      
      kernel_input <- total_input[,group == fold][group != fold,]
      
      f_i_hats <- colSums(kernel(kernel_input))
      
      fold_l_cv[fold] <- sum(log(1/(h * N_i) * f_i_hats))
    }
    
    cv_sum[h_index] <- sum(fold_l_cv)
    
  }
  
  cv_result <- data.frame(h = h_grid, cv_result = cv_sum)
  best_h <- h_grid[which.max(cv_sum)]
  
  return(list(cv_result = cv_result, best_h = best_h))
}



cv_currying <- function(x, h_grid, folds, kernel = epanechnikov, seed = NA) {
  
  N <- length(x)
  cv_sum <- numeric(length(h_grid))
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  group <- sample(rep(1:folds, length.out = N))
  
  for (h_index in seq_along(h_grid)){
    
    h <- h_grid[h_index]
    fold_l_cv <- numeric(folds)
    
    
    for (fold in seq_along(1:folds)){
      other_sets <- x[group != fold]
      hold_out_set <- x[group == fold]
      
      outer(other_sets, -hold_out_set/h, '%*%')
      
      f_i_hats <- numeric(length(hold_out_set))
      N_i <- length(other_sets)
      
      epa_pre_calc <- function(other_sets, h){
        return(
          function(hold_out_set){
            return(epanechnikov((hold_out_set - other_sets) / h))
          })
      }
      
      epa_calc <- epa_pre_calc(other_sets, h)
      
      fold_l_cv[fold] <- sum(log(1/(h * N_i) * colSums(sapply(hold_out_set, epa_calc))))
    }
    
    cv_sum[h_index] <- sum(fold_l_cv)
    
  }
  
  cv_result <- data.frame(h = h_grid, cv_result = cv_sum)
  best_h <- h_grid[which.max(cv_sum)]
  
  return(list(cv_result = cv_result, best_h = best_h))
}





testCV <- function(x, h_grid, folds, kernel = epanechnikov, B, seed = NA){
  
  n <- length(x)
  PEcv <- vector("list", B)
  tmp <- numeric(n)
  
  for (b in 1:B){
    group <- sample(rep(1:folds, length.out = n))
    
    for (h_index in seq_along(h_grid)){
      h <- h_grid[h_index]
      fold_l_cv <- numeric(folds)
      
      for (i in 1:folds){
        # Define vector to store f_i_hat for each x_i in hold_out_set
        f_i_hats <- numeric(length(x[group == fold]))
        N_i <- length(x[group != fold])
        
        for (x_i in seq_along(x[group == fold])){
          f_i_hat <- 1/(h * N_i) * sum(epanechnikov((x[group == fold][x_i] - x[group != fold]) / h))
          f_i_hats[x_i] <- f_i_hat
        }
        fold_l_cv[fold] <- sum(log(f_i_hats))
      }
    }
  }
}


#------------ Testing ------------
# x <- rnorm(500)
# #cv_naive(x, h_grid = seq(0.01, 1.5, by = 0.01), folds = 15, kernel = epanechnikov, seed = 123)
# cv_outer(x, h_grid = seq(0.01, 1.5, by = 0.01), folds = 15, kernel = epanechnikov, seed = 123)
# #cv_Outer(x, h_grid = seq(0.01, 1.5, by = 0.01), folds = 15, kernel = epanechnikov, seed = 123)
# 
# #profvis::profvis(cv_outer(x, h_grid = seq(0.01, 1, by = 0.05), folds = 10, kernel = epanechnikov, seed = 123))
# 
# # bench::mark(cv_outer(x, h_grid = seq(0.01, 1.5, by = 0.01), folds = 15, kernel = epanechnikov, seed = 123), 
# #             cv_Outer(x, h_grid = seq(0.01, 1.5, by = 0.01), folds = 15, kernel = epanechnikov, seed = 123),
# #             check = F)







