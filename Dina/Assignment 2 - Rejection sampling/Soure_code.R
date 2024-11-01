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
library(foreach)
library(doParallel)

## Getting the data ###########################################
poisson_data <- read_csv("C:/Users/birgi/Documents/comp_stat/Dina/Assignment 2 - Rejection sampling/poisson.csv")


## Density function and zx term ###################
zx <- sum(poisson_data$x*poisson_data$z)

tar_dens <- function(y){
  sapply(y, function(y) exp(y * zx - sum(exp(y*poisson_data$x))))
}


## Function for optimal mu and sigma for gauss envelope ###################
optim_s_alpha <- function(mu){
  y_seq <- seq(0,1, 0.001)
  sigma <- function(s){
    alpha_star <- min(dnorm(y_seq, mu_opt, s)/tar_dens(y_seq))
    return(-alpha_star)
  }
  
    sigma_opt <- optimize(sigma, c(0,1))$minimum
    alpha_s <- -optimize(sigma, c(0,1))$objective
    
    return(list(sigma = sigma_opt, alpha = alpha_s))
  }
 
## log of target density ###################

tar_dens_log <- function(y){
  sapply(y, function(y) log(tar_dens(y)))
}

## Gaussian rejection sampler - vectorized  ###################

gauss_rej <- function(N, mu, sigma, alpha, dens) {
      y0 <- rnorm(N, mean = mu, sd = sigma)
      u <- runif(N)
      reject <- u > 1/(sqrt(2*pi*sigma)) * alpha * dens(y0) / 
        dnorm(y0, mean = mu, sd = sigma)
      
      
     y <- y0[!reject]

  return(y)
}

## Function factory for getting fast sampler ############################

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

## Fast samplers ############################

gaus_rej_f <- new_rejection_sampler(gauss_rej)



#The derivative of log(tar_dens) in point y
tar_dens_log_difference <- function(y){
  diff <- sapply(y, function(y) zx - sum(poisson_data$x * exp(y * poisson_data$x)))
  diff
}

#Calculating H_i's
H_i <- function(as, bs, zs, n) {
  1 / as * exp(bs) * (exp(as * zs[2:(n+1)]) - exp(as * zs[1:n]))
}

z_fun <- function(a, b, start_point, end_point){
  # Defining z:
  z <- numeric(length(a) + 1)  # length(a) + 1 because we have z_0 and z_m
  
  # Set the boundary points
  z[1] <- start_point  # z_0
  z[length(z)] <- end_point  # z_m
  
  # Calculating intermediate z_i values
  for (i in 1:(length(a) - 1)) {
    z[i + 1] <- (b[i + 1] - b[i]) / (a[i] - a[i + 1])
  }
  
  # Return the resulting z vector
  return(z)
}


cor_interval <- function(points, interval_sum){
  # Find which interval each point belongs to
  intervals <- findInterval(points, interval_sum, left.open = T)
  
  # Create a matrix of 0's with rows = length(points) and columns = length(cumsum_vec)
  result_matrix <- matrix(0, nrow = length(points), ncol = length(interval_sum))
  
  # Set the appropriate interval to 1 for each point
  result_matrix[cbind(seq_along(points), intervals + 1)] <- 1
  
  return(result_matrix)
}


extract_nonzero_indices <- function(mat) {
  apply(mat, 1, function(row) which(row != 0))
}

adap_log_rej_sampler <- function(N, breakpoints, log_tar_dens, log_tar_dens_diff){
  
  #Calculating a's, b's and z's given the breakpoints:
  a <- log_tar_dens_diff(breakpoints)
  #Checking that a is non-zero, and that they are unique
  if( anyDuplicated(a) > 0 || any(a == 0) ){
    stop("\n
    Error: The a-coefficients are required to be non-zero and unique. 
    Change the choice of breakpoints."
    )
  }
  b <- log_tar_dens(breakpoints) - a * breakpoints
  z <- z_fun(a, b, 0, 5)
  
  #Now calculating R_i and Q_i's
  R_i <- H_i(a, b, z, length(a))
  Q_i <- cumsum(R_i)
  
  #Defining d:
  d <- Q_i[length(Q_i)]
  
  #Calculating dq, and checking whether it falls into the right intervals:
  dq <- d * runif(N)
  I_matrix <- cor_interval(dq, Q_i)
  
  #Now we transform this matrix into the one containing dq, and extract the indices of the rows, so we know which a_i, b_i etc. 
  #correspond to each dq:
  dq_matrix <- I_matrix * dq
  interval_indicator <- extract_nonzero_indices(I_matrix)
  
  
  #Now we calculate the necessary Q_i-1's, a_i's etc.
  a_i <- a[interval_indicator]
  b_i <- b[interval_indicator]
  Q_i_minus <- c(0, Q_i[-length(Q_i)])[interval_indicator]
  z_i_minus <- z[-length(z)][interval_indicator]
  
  # Finally we sample from the proposal distribution
  x <- 1 / a_i * log(a_i * exp(-b_i) * (dq - Q_i_minus) + exp(a_i*z_i_minus))
  
  accept <- logical(N)
  u <- runif(N)
  
  for (i in 1:length(a)){ #1 to m
    group_i <- as.logical(I_matrix[, i]) # Index to select observations belonging to interval i
    
    accept[group_i] <- u[group_i] <= exp(log_tar_dens(x[group_i]) - a[i] * x[group_i] - b[i])
  }
  return(x[accept])
  #return(list(x[accept], a = a, b = b, z = z))
}






## For plotting envelopes ############################

affine_var <- function(breakpoints){
  a <- tar_dens_log_difference(breakpoints)
  b <-  tar_dens_log(breakpoints) - a * breakpoints
  z <- z_fun(a, b, 0, 5)
  
  return(list(a = a, b = b, z = z))
}

exp_affine_values <- function(y, a, b) {
  # Create a matrix to hold the affine values
  val_matrix <- matrix(NA, nrow = length(a), ncol = length(y))
  
  # Calculate the affine values for each a and b
  for (i in 1:length(a)) {
    val_matrix[i, ] <- exp(a[i] * y + b[i])  # Vectorized operation for each row
  }
  
  return(val_matrix)
}

## Other version of the envelope  ############################

####################################################
########### Piecewise linear Envelope ##############
####################################################

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

piece_lin_rejec_samp <- function(N, ys) {
  
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
  for(i in 1:N){
    # Finding the interval x_i belongs to
    int <- which(I[i,] == 1)
    
    # Taking the inverse cdf
    x[i] <- log((u0[i] - Q[int]) * as[int] * exp(- bs[int]) + exp(as[int] * zs[int])) / as[int]
    
    # Acceptance step
    accept[i] <- u[i] <=  tar_dens(x[i]) / exp(as[int] * x[i] + bs[int])
  }
  
  return(x[accept])
}

aff_rej_n <- new_rejection_sampler(piece_lin_rejec_samp)


## Slow envelope: ###########################################
piece_lin_rejec_samp_s <- function(N, ys) {
  
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
  I <- matrix(FALSE, nrow = N, ncol = n)
  
  for (i in 1:N) {
    for (j in 1:n) {
      I[i, j] <- (u0[i] > Q[j]) && (u0[i] <= Q[j + 1])
    }
  }
  
  x <- numeric(N)
  accept <- logical(N)
  
  for (i in 1:N) {
    # Finding the interval x_i belongs to
    int <- which(I[i,] == TRUE)
    
    # Taking the inverse cdf
    if (length(int) > 0) {  # Check if the interval is found
      x[i] <- log((u0[i] - Q[int]) * as[int] * exp(-bs[int]) + exp(as[int] * zs[int])) / as[int]
      
      # Acceptance step
      accept[i] <- u[i] <= tar_dens(x[i]) / exp(as[int] * x[i] + bs[int])
    } else {
      x[i] <- NA  # Set to NA if no interval found
      accept[i] <- FALSE
    }
  }
  
  return(x[accept])
}

aff_rej_s <- new_rejection_sampler(piece_lin_rejec_samp_s)
aff_rej_n(10, seq(0,0.35, length.out = 10))
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
  results <- foreach(i = 1:N, .combine = rbind, .export = c("tar_dens", "zx", "poisson_data")) %dopar% {
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
  
  
  return(x[accept])
}

aff_rej_par <- new_rejection_sampler(piece_lin_rejec_samp_par)


## Rcpp version of tar_dens ############################

sourceCpp("~/comp_stat/Dina/Assignment 2 - Rejection sampling/Rcpp_loop.cpp")

tar_dens(seq(0, 0.35, length.out = 10))
tar_den_cpp(seq(0, 0.35, length.out = 10), )

tar_dens_cpp(seq(0, 0.35, length.out = 10), poisson_data$x, zx)
tar_dens_log_difference <- function(y){
  diff <- sapply(y, function(y) zx - sum(poisson_data$x * exp(y * poisson_data$x)))
  diff
}

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

piece_lin_rejec_samp <- function(N, ys) {
  
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
  for(i in 1:N){
    # Finding the interval x_i belongs to
    int <- which(I[i,] == 1)
    
    # Taking the inverse cdf
    x[i] <- log((u0[i] - Q[int]) * as[int] * exp(- bs[int]) + exp(as[int] * zs[int])) / as[int]
    
    # Acceptance step
    accept[i] <- u[i] <=  tar_dens(x[i]) / exp(as[int] * x[i] + bs[int])
  }
}


