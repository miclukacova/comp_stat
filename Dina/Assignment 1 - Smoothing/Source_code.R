##################################################
################### Packages #####################
##################################################

library(bench)
library(ggplot2)
library(tidyverse)
library(testthat)
library(profvis)
library(gridExtra)
theme_set(theme_bw())
library(CSwR)

## Function to obtain pilot bw:###################
r_hat <- function(x){
  r_hat <- 2.34*min(sd(x), IQR(x)/1.34)*length(x)^(-1/5)
  return(r_hat)
}

## AMISE Optimal bw-function ###################
h_N_am <- function(x){
  #Defining parameters:
  r <- r_hat(x)
  n <- length(x)
  
  #Calculating the f-tilde L2 norm
  f_norm <- 9/4 * 1/(r^6 * n^2) * sum(pmax(0, outer(x, x, pmin) - outer(x, x, pmax) + 2 * r))
  
  #Calculating the optimal bandqwidth using the f_norm:
  h_N <- ( (3/5) / (f_norm * 1/25) )^(1/5) * n^(-1/5)
  
  return(h_N)
}

## Naive CV bw-function ###################
cv_bw_M <- function(x, k = 10, h = seq(0.03, 10, by = 0.01)){
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
      #if(any(f_hat_i < 0 )) {print(f_hat_i)} #Check for negative values of f_hat
    }
    
    #Calculating log likelihood
    bw_l[j] <- sum(log(f_hat_i))
    #print(c(bw_l[j], h[j]))
    
  }
  
  #returning max likelihood:
  return(h[which.max(bw_l)])
}

## Vectorized CV bw-function ###################
bw_cv_vec <- function(x, k = 10, h = seq(0.03, 10, by = 0.01)) {
  # Creating the kernel function:
  kernel <- function(z) (abs(z) <= 1) * 3 / 4 * (1 - z^2)
  # Creating the h's to loop over:
  h <- h
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
      
      #if(any(f_hat_i < 0 )) {print(f_hat_i)} #Check for negative values of f_hat
      
    }
    # Calculate log-likelihood for current h
    ll[j] <- sum(log(f_hat))
    #print(c(ll[j], h[j]))
  }
  return(h[which.max(ll)])
}
## kernel function ###########################################

kernel <- function(z) (abs(z) <= 1) * 3 / 4 * (1 - z^2) 


## no of observations in each bin ############################
kern_bin <- function(x, l, u, B) {
  w <- numeric(B)
  delta <- (u - l) / (B - 1)
  for (j in seq_along(x)) {
    i <- floor((x[j] - l) / delta + 0.5) + 1
    w[i] <- w[i] + 1
  }
  w / sum(w)
}


## Kernel density estimator ############################

kern_dens_bin <- function(x, h = 1, m = 512, B, normalize = TRUE, k = 10) {
  #Selecting the right bw:
  if(h == "AMISE"){
    h <- h_N_am(x)
  } else if (h == "CV"){
    h <- bw_cv_vec(x, k)
  } else if (is.numeric(h) == FALSE){
    stop("Invalid bandwidth. Please choose a double or 'AMISE' or 'CV'")
  }
  
  #finding the range of the data, defining our gridpoints and y
  rg <- range(x)
  y <- numeric(m)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  #Normalizing if relevant:
  if(normalize){
    h_n <- h * sqrt(5) #Variance of Epanechnikov kernel
  }
  h_n <- h
  #Calculating delta and finding the center points.
  delta <- (rg[2] - rg[1]) / (B - 1)
  cent_points <- seq((rg[1] + delta/2), (rg[2] + delta/2), len = B) 
  #Envoking the kern_bin function to find the weights/n_j's
  nj <- kern_bin(x, rg[1], rg[2], B)
  
  for (i in seq_along(xx)) {
    y[i] <- mean(nj %*% (kernel((xx[i] - cent_points)/h_n)))/(h_n)
  }
  list(x = xx, y = y, h = h)
}



###### S3 Classes ###########################################

#-------------------------------------------------------------------------------
#                                   kernel density class                                  #
#-------------------------------------------------------------------------------

#  class
kernel_density <- function(x, h, B, m = 512, normalize = TRUE, k = 10) {
  structure(
    list(
      density = kern_dens_bin(x = x, h = h, m = m, B = B, normalize = normalize),
      data = x),
    class = "kernel_density"
  )
}


# Summary method
summary.kernel_density <- function(object) {
  return(head(data.frame(x = object$density$x, y = object$density$y)))
}


# Print method for objects of class "kernel_density"
print.kernel_density <- function(object) {
  cat("Kernel density\n")
  cat("Bandwidth h:\n")
  print(object$density$h)
  
  cat("\nDensity estimates:\n")
  
  # Print the first 5 x and y values individually without converting to data.frame
  x_values <- object$density$x[1:5]
  y_values <- object$density$y[1:5]
  
  # Display x and y values together in a clear format
  for (i in seq_along(x_values)) {
    cat(sprintf("x: %f, y: %f\n", x_values[i], y_values[i]))
  }
}


# Plot method for objects of class "kernel_density"
plot.kernel_density <- function(object) {
  # Check if the object is of the correct class
  if (!inherits(object, "kernel_density")) {
    stop("Object must be of class 'kernel_density'")
  }
  
  # Extract the density values from the object
  e_df <- data.frame(x = object$density$x, y = object$density$y)
  
  # Plot using ggplot2
  ggplot(e_df, aes(x = x)) +
    geom_line(aes(y = y), color = "#4f7942", size = 1) +
    geom_histogram(data = data.frame(x = object$data), aes(x = x, y = ..density..), 
                   fill = "coral", color = "coral1", alpha = 0.2, bins = 40) +
    labs(x = "Input values - x", y = "Density") +
    theme_bw()
}


plot_data <- function(x) {
  UseMethod("plot_data")
}


plot_data.kernel_density <- function(object) {
  x <- object$density$x
  y <- object$density$y
  
  kd_plot <- data.frame(x = x, y = y)
  
  return(kd_plot)
}



## Paralellizing the CV bw-function ###################

bw_cv_vec_parallel <- function(x, k = 10, h = seq(0.03, 10, by = 0.01)) {
  # Creating the kernel function:
  kernel <- function(z) (abs(z) <= 1) * 3 / 4 * (1 - z^2)
  
  # Defining n, the length of the dataset:
  n <- length(x)
  
  # Randomly dividing the number of indices so that we get k-groups
  group <- sample(rep(1:k, length.out = n))
  
  # Setting up parallel backend
  cl <- makeCluster(detectCores() - 1)  # Leave one core free
  registerDoParallel(cl)
  
  # Use foreach to parallelize the computation of log-likelihoods for each h
  ll <- foreach(j = 1:length(h), .combine = c) %dopar% {
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
    sum(log(f_hat))
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Find the index of the maximum log-likelihood
  return(h[which.max(ll)])
}


## Trying to implement the inner for loop in C++ ###################
library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/comp_stat/Dina/Assignment 1 - Smoothing/CV_bw_rcpp.cpp")


bw_cv_rcpp <- function(x, k = 10, h = seq(0.03, 10, by = 0.01)) {
  # Defining n, the length of the dataset:
  n <- length(x)
  # Randomly dividing the number of indices so that we get k groups
  group <- sample(rep(1:k, length.out = n))
  
  # For each h, calculate the log-likelihood estimate
  ll <- calculate_log_likelihood(x, group, h, k)
  
  # Find the h that maximizes the log-likelihood
  optimal_h <- h[which.max(ll)]
  
  return(optimal_h)
}





