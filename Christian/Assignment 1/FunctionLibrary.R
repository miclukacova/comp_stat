###############################################
#                                             #
#              Function library               #
#                                             #
###############################################



#-----------------Function to calculate the Epanechnikov kernel-----------------
epanechnikov <- function(x) {
  (abs(x) <= 1) * 3/4 * (1 - x^2)
}
#-------------------------------------------------------------------------------




#------------------------Function to prepare the kernel-------------------------
kernel_prep <- function(x, h, n){
  if (h <= 0) stop("Bandwidth must be positive")
  
  # Get the range and length of the input vector
  rg <- range(x)
  N <- length(x)
  
  # Define the grid points corresponding to the gridpoints from density()
  gridpoints <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = n)
  y <- numeric(n)
  
  return(list(range = rg, gridpoints = gridpoints, y = y, N = N))
}
#-------------------------------------------------------------------------------





#---------------------------Function to determine bins--------------------------
kern_bin <- function(x, l, u, B) {
  
  # Create vectors to store weights and centers
  w <- numeric(B)
  
  # Define interval length
  delta <- (u - l) / (B - 1)
  
  # Create vector of centers
  centers <- seq(l + delta/2, u - delta/2, length.out = B)
  
  # Loop through each data point
  for (j in seq_along(x)) {
    i <- floor((x[j] - l) / delta + 0.5) + 1
    w[i] <- w[i] + 1
  }
  return(list(weights = w, centers = centers))
}



kern_bin_Rcpp <- function(x, l, u, B) {
  
  # Define interval length
  delta <- (u - l) / (B - 1)
  
  # Create vector of centers
  centers <- seq(l + delta/2, u - delta/2, length.out = B)
  
  # Loop through each data point
  w <- kernbinloopRcpp(B = B, x = x, l = l, delta = delta)
  
  return(list(weights = w, centers = centers))
}


#-------------------------------------------------------------------------------


calc_epa <- function(weights, centers, h){
  return(function(gridpoints){
    sum(weights * epanechnikov((gridpoints - centers) / h))
  })
} 
  








#-------------------------------Calculate sigma_tilde---------------------------
sigma_tilde <- function(x){
  min(sd(x), IQR(x)/1.34)
}
#-------------------------------------------------------------------------------




#-----------------------Calculate Epanechnikov pilot bandwidth------------------
epanechnikov_pilot_bandwidth <- function(x) {
  return((40 * sqrt(pi))^(1/5) * sigma_tilde(x) * length(x)^(-1/5))
}
#-------------------------------------------------------------------------------




#-----------Calculate difference between min and max in x pilot bandwidth-------
min_max_diff_sum <- function(x, r_hat){
  min_matrix <- outer(x,x, pmin)
  max_matrix <- outer(x,x, pmax)
  
  return(sum(pmax(0, 2*r_hat + min_matrix - max_matrix)))
}
#-------------------------------------------------------------------------------


#--------Calculate difference between min and max in x pilot bandwidth fast-----
min_max_diff_sum_fast <- function(x, r_hat){
  abs_val_matrix <- - abs(outer(x,x, "-"))
  
  return(sum(pmax(0,2*r_hat + abs_val_matrix)))
}
#-------------------------------------------------------------------------------


#-------Calculate difference between min and max in x pilot bandwidth faster----
min_max_diff_sum_faster <- function(x, r_hat){
  abs_val_matrix <- - abs(Outer(x,x, "-"))
  
  return(sum(pmax(0,2*r_hat + abs_val_matrix)))
}


min_max_diff_sum_Rcpp <- function(x, r_hat){
  abs_val_matrix <- - abs(Outer(x,x, "-"))
  
  return(matrix_sum_pmaxRcpp(r_hat, abs_val_matrix))
}
#-------------------------------------------------------------------------------





#-----------------Calculate f_tilde_norm2 for Epanechnikov kernel---------------
epanechnikov_f_tilde_norm2 <- function(x){
  
  n <- length(x)
  
  r_hat <- epanechnikov_pilot_bandwidth(x)
  
  double_sum <- min_max_diff_sum_faster(x, r_hat)
  
  return(9/4 * 1/(n^2 * r_hat^6) * double_sum)
  
}

epanechnikov_f_tilde_norm2_rcpp <- function(x){
  
  n <- length(x)
  
  r_hat <- epanechnikov_pilot_bandwidth(x)
  
  double_sum <- min_max_diff_sum_Rcpp(x, r_hat)
  
  return(9/4 * 1/(n^2 * r_hat^6) * double_sum)
  
}
#-------------------------------------------------------------------------------





