#-------------------------------------------------------------------------------
# Binned kernel estimator
#-------------------------------------------------------------------------------

kern_dens_bin <- function(x, h, m = 512, B = 100, norm = TRUE, bin = kern_bin) {
  rg <- range(x)
  delta <- (rg[2] - rg[1]) / (B - 1)
  c <- seq(rg[1] + delta/2, rg[2] - delta/2, length.out = B)
  n <- bin(x, rg[1], rg[2], B)
  
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