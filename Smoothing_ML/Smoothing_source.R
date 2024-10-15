library(bench)
library(ggplot2)
library(tidyverse)
library(testthat)
library(bench)
library(profvis)
theme_set(theme_bw())

########### Functions ##############

#-------------------------------------------------------------------------------
# Naive density estimation
#-------------------------------------------------------------------------------

kern_dens1 <- function(x, h, m = 512, norm = TRUE) {
  
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m)
  
  # we normalize now, so that the gridpoints are equal to r's density function
  if(norm){
    h <- sqrt(5) * h
  }
  
  for (i in seq_along(xx)) {
    condition <- abs(xx[i] - x) <= h
    y[i] <- sum((1 - (xx[i] - x)^2 / h^2) * condition)
  }
  y <- 3 * y / (4 * h * length(x))
  list(x = xx, y = y)
  
}

#-------------------------------------------------------------------------------
# Function computing the percentage of observations in each bin
#-------------------------------------------------------------------------------

kern_bin <- function(x, l, u, B) {
  w <- numeric(B)
  delta <- (u - l) / (B - 1)
  for (j in seq_along(x)) {
    i <- floor((x[j] - l) / delta + 0.5) + 1
    w[i] <- w[i] + 1
  }
  w / sum(w)
}

#-------------------------------------------------------------------------------
# Binned kernel estimator
#-------------------------------------------------------------------------------

kern_dens_bin <- function(x, h, m = 512, B = 100, norm = TRUE) {
  rg <- range(x)
  delta <- (rg[2] - rg[1]) / (B - 1)
  c <- seq(rg[1] + delta/2, rg[2] - delta/2, length.out = B)
  n <- kern_bin(x, rg[1], rg[2], B)
  
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m)
  
  # we normalize now, so that the gridpoints are equal to r's density function
  if(norm){
    h <- sqrt(5) * h
  }
  
  for (i in seq_along(xx)) {
    condition <- abs(xx[i] - c)  <=  h
    y[i] <- sum(n * (1 - (xx[i] - c)^2 / h^2) * condition)
  }
  
  y <- 3 * y / (4 * h)
  list(x = xx, y = y)
}

#-------------------------------------------------------------------------------
# AMISE based optimal bandwidth finder
#-------------------------------------------------------------------------------

AMISE_bw <- function(x) {
  
  n <- length(x) 
  
  # Silverman bandwidth
  norm_h <- (40 * sqrt(pi))^(1/5) * min(sd(x), IQR(x) / (1.34)) * n^(-1/5)
  
  # Pilot density estimate
  norm_p_f <- 0
  for(i in seq_along(x)){
    int <- pmax(0, 2 * norm_h - abs(x - x[i]))  
    norm_p_f <- norm_p_f + sum(int)
  }
  norm_p_f <- 9 * norm_p_f /(4*n^2 * norm_h^6)
  
  # AMISE bandwidth
  (15 / norm_p_f)^(1/5) * n^(-1/5)
}

#-------------------------------------------------------------------------------
# CV based optimal band width finder
#-------------------------------------------------------------------------------

cv_bw_l <- function(x, k = 5, h = seq(0.01, 2, by = 0.01)) {
  
  # Partition
  n <- length(x)
  groups <- sample(rep_len(seq(1,k), length.out = n), replace = FALSE)
  
  # number of obs in each group
  N <- sapply(seq(1,k), function(x) sum(groups == x), simplify = TRUE)
  
  # bandwidth log likelihoods 
  bw_l <- numeric(length(h))
  
  # the kernel density estimates
  f_hat_i <- numeric(length(x))
  
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

#-------------------------------------------------------------------------------
# MSE of density
#-------------------------------------------------------------------------------

mse_dens <- function(est, true) {
  mean((est - true)^2)
}



########### Class ##############

density_object <- function(x, h = NULL) {
  structure(
    list(
      x = x,
      h = h,
      cv_bw = cv_bw_l_fast(x),
      amise_bw = AMISE_bw(x)
    ),
    class = "density_object"
  )
}

kern_dens <- function(x) {
  UseMethod("kern_dens")
}

kern_dens <- function(object, h = "CV", binned = FALSE) {
  if(is.null(h)){
    h <- object$cv_bw
  }
  else if(h == "CV"){
    h <- object$cv_bw
  }
  else if(h == "AMISE"){
    h <- object$amise_bw
  }
  if(binned){
      kern_dens_bin(object$x, h)
    }
  kern_dens1(object$x, h)
}
