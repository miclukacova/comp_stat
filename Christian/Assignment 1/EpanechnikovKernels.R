###############################################
#                                             #
#     Epanechnikov kernel implementations     #
#                                             #
###############################################


#-------------------------------Loop based kernel-------------------------------
kern_dens_loop <- function(x, h, n = 512, normalization = TRUE) {
  
  if (h <= 0) stop("Bandwidth must be positive")
  
  # Prepare the kernel
  prep <- kernel_prep(x, h, n)
  gridpoints <- prep$gridpoints
  y <- prep$y
  N <- prep$N
  
  # Check if the normalization is requested
  if (normalization){
    # Multiply bandwidth by sqrt() inverse kernel variance if normalization is requested
    h <- h * sqrt(5)
  }
  
  const <- h * N
  
  # The computation is done using nested for-loops. The outer loop
  # is over the grid points, and the inner loop is over the data points.
  for (i in seq_along(gridpoints)) {
    for (j in seq_along(x)) {
      y[i] <- y[i] + epanechnikov((gridpoints[i] - x[j]) / h)
    }
  }
  
  # Normalize the density estimate
  y <- y / const
  
  list(x = gridpoints, y = y)
}
#-------------------------------------------------------------------------------




#-------------------------------Vectorized kernel-------------------------------
kern_dens_vec <- function(x, h, n = 512, normalization = TRUE) {
  
  if (h <= 0) stop("Bandwidth must be positive")
  
  # Prepare the kernel
  prep <- kernel_prep(x, h, n)
  gridpoints <- prep$gridpoints
  y <- prep$y
  N <- prep$N
  
  # Check if the normalization is requested
  if (normalization){
    # Multiply bandwidth by sqrt() inverse kernel variance if normalization is requested
    h <- h * sqrt(5)
  }
  
  const <- h * N
  
  # The inner loop from 'kern_dens_loop()' has been vectorized, and 
  # only the outer loop over the grid points remains. 
  
  for (i in seq_along(gridpoints)) {
    y[i] <- sum(epanechnikov((gridpoints[i] - x) / h)) / const
  }
  
  list(x = gridpoints, y = y)
}
#-------------------------------------------------------------------------------




#---------------------------------Binned kernel---------------------------------

kern_dens_bin <- function(x, h, n = 512, B = 100, normalization = TRUE, fast = TRUE) {
  
  if (h <= 0) stop("Bandwidth must be positive")
  
  # Prepare the kernel
  prep <- kernel_prep(x, h, n)
  gridpoints <- prep$gridpoints
  y <- prep$y
  N <- prep$N
  range <- prep$range
  
  # Check if the normalization is requested
  if (normalization){
    # Multiply bandwidth by sqrt() inverse kernel variance if normalization is requested
    h <- h * sqrt(5)
  }
  
  const <- h * N
  
  # Define the binning characteristics
  bins <- kern_bin(x = x, l = range[1], u = range[2], B = B)
  weights <- bins$weights
  centers <- bins$centers
  
  # Calculate binned kernel density estimates
  if (fast){
    y <- outer(gridpoints, centers, FUN = function(x,y) epanechnikov((x - y) / h))
    y <- c(y %*% weights)
  } else {
    for (i in seq_along(gridpoints)) {
      y[i] <- sum(weights * epanechnikov((gridpoints[i] - centers) / h))
    }
  }
  
  y <- y / const
  
  
  list(x = gridpoints, y = y)
}





kern_dens_bin_Rcpp <- function(x, h, n = 512, B = 100, normalization = TRUE, fast = TRUE) {
  
  if (h <= 0) stop("Bandwidth must be positive")
  
  # Prepare the kernel
  prep <- kernel_prep(x, h, n)
  gridpoints <- prep$gridpoints
  y <- prep$y
  N <- prep$N
  range <- prep$range
  
  # Check if the normalization is requested
  if (normalization){
    # Multiply bandwidth by sqrt() inverse kernel variance if normalization is requested
    h <- h * sqrt(5)
  }
  
  const <- h * N
  
  # Define the binning characteristics
  bins <- kern_bin_Rcpp(x = x, l = range[1], u = range[2], B = B)
  weights <- bins$weights
  centers <- bins$centers
  
  # Calculate binned kernel density estimates
  if (fast){
    y <- c(outer(gridpoints, centers, FUN = function(x,y) epanechnikov((x - y) / h)) %*% weights)
  } else {
    for (i in seq_along(gridpoints)) {
      y[i] <- sum(weights * epanechnikov((gridpoints[i] - centers) / h))
    }
  }
  
  y <- y / const
  
  
  list(x = gridpoints, y = y)
}

#-------------------------------------------------------------------------------





# Use sparcity?









