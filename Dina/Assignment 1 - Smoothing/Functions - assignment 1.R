cv_bw_M <- function(x, k = 5, h = seq(0.01, 2, by = 0.01)){
  #Partition
  n <-length(x)
  groups <- sample(rep(1:k, length.out = n))
  
  #Number of observations in each group:
  N <- sapply(seq(1,k), function(x) sum (groups == x), simplify = T)
  
  #bandwidth log likelihoods:
  bw_l <- numeric(length(h))
  
  for (j in seq_along(h)) {
    #kernel density estimates:
    f_hat_i <- numeric(length(x))
    
    #Calculating CV density estimates
    for (i in seq_along(x)) {
      k <- groups[i]
      x_m_i <- x[groups != k]
      condition <- abs((x[i] - x_m_i)/h[j]) <= 1
      f_hat_i[i] <- 3 / (N[k] * 4 * h[j]) * sum((1-(x[i] - x_m_i)^2/h[j]^2) * condition)
      
    }
    
    #Calculating log likelihood
    bw_l[j] <- sum(log(f_hat_i))
    #print(c(bw_l[j], h[j]))
    
  }
  
  #returning max likelihood:
  return(h[which.max(bw_l)])
}


bw_cv_vec <- function(x, k = 5, h = seq(0.01, 2, by = 0.01)) {
  # Creating the kernel function:
  kernel <- function(z) (abs(z) <= 1) * 3 / 4 * (1 - z^2)
  
  # Creating the h's to loop over:
  h <- seq(0.01, 2, 0.02)
  
  # Defining n, the length of the dataset:
  n <- length(x)
  
  # Randomly dividing the number of indices so that we get k-groups
  group <- sample(rep(1:k, length.out = n))
  
  # Creating our return vectors:
  ll <- numeric(length(h))  # To store log-likelihood for each h
  
  # For each h, calculate the log-likelihood estimate:
  for (j in 1:length(h)) {
    f_hat <- numeric()  # Reset f_hat for each h
    
    # For each fold, calculate the f_hat_i's
    for (i in 1:k) {
      # Number of indices not in group i
      N_i <- length(x[group != i])
      
      # Creating a matrix consisting of all points (x_i - x_j)/h
      x_matrix <- 1 / h[j] * outer(x[group != i], x[group == i], "-")
      
      # Evaluating the kernel in each point
      kerns <- kernel(x_matrix)
      
      # Summing the kernel values over columns
      kern_sums <- colSums(kerns)
      
      # Calculating f_hat for the given group
      f_hat_i <- 1 / (h[j] * N_i) * kern_sums
      
      # Store the results
      f_hat <- c(f_hat, f_hat_i)
      
    }
    
    # Calculate log-likelihood for current h
    ll[j] <- sum(log(f_hat))
    #print(c(ll[j], h[j]))
  }
  
  # Return the h that minimizes the log-likelihood
  return(h[which.max(ll)])
}









