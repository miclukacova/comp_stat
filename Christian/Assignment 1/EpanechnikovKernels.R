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

kern_dens_bin <- function(x, h, n = 512, B = 100, normalization = TRUE) {
  
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
  # Calculate binned kernel density estimates
  
  y <- c(outer(gridpoints, centers, FUN = function(x,y) epanechnikov((x - y) / h)) %*% weights)
  
  # for (i in seq_along(gridpoints)) {
  #     y[i] <- sum(weights * epanechnikov((gridpoints[i] - centers) / h))
  # }
  
  y <- y / const
  
  
  list(x = gridpoints, y = y)
}





kern_dens_bin_Rcpp <- function(x, h, n = 512, B = 100, normalization = TRUE) {
  
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
  y <- c(outer(gridpoints, centers, FUN = function(x,y) epanechnikov((x - y) / h)) %*% weights)
  
  # for (i in seq_along(gridpoints)) {
  #     y[i] <- sum(weights * epanechnikov((gridpoints[i] - centers) / h))
  # }
  
  
  y <- y / const
  
  
  list(x = gridpoints, y = y)
}

#-------------------------------------------------------------------------------





# Object

epanechnikovKernel <- function(x, 
                               h = "Silverman", 
                               n = 512, 
                               normalization = TRUE,
                               h_grid = seq(0, 2, 0.05),
                               kernel = "binned",
                               seed = NA, 
                               B = 100,
                               folds = 10,
                               reps = 50) {
  
  if (h <= 0) stop("Bandwidth must be positive")
  
  bandwidth_method <- h
  
  if (h == "AMISE"){
    h <- AMISE_epa_bandwidth_rcpp(x)
  } else if (h == "CV") {
    h <- CV(x = x, h_grid = h_grid, folds = folds, seed = seed, reps = reps)
  } else if (h == "Silverman"){
    h <- silverman(x)
  }
  
  if (kernel == "vectorized"){
    kernel_est <- kern_dens_vec(x = x, h = h, n = n, normalization = normalization)
  } else if (kernel == "binned") {
    kernel_est <- kern_dens_bin_Rcpp(x = x, h = h, n = n, B = B, normalization = normalization)
  }
  
  output <- structure(
    list(
      x = kernel_est$x,
      y = kernel_est$y,
      h = h,
      bandwidth_method = bandwidth_method,
      kernel_method = kernel,
      data = x
    ),
    class = "epanechnikovKernel"
  )
}

plot.epanechnikovKernel <- function(object, 
                                    bins = 30, 
                                    color = "white", 
                                    fill = "steelblue", linecolor = "red", lty = 1, lwd = 0.5, ...) {
  
  data_df <- data.frame(object$data)
  kernel_est_df <- data.frame(x = object$x, y = object$y)
  
  ggplot(data = data_df, aes(x = x)) + geom_histogram(aes(y = ..density..), color = color, fill = fill) +
    geom_line(data = kernel_est_df, aes(x = x, y = y), color = linecolor, lty = lty, lwd = lwd) +
    ggtitle("Epanechnikov Kernel Density Estimate")
}

print.epanechnikovKernel <- function(object, ...) {
  cat("Epanechnikov Kernel Density Estimate\n")
  cat("Bandwidth: ", object$h, "\n")
  cat("Bandwidth method: ", object$bandwidth_method, "\n")
  cat("Kernel method: ", object$kernel_method, "\n")
}



