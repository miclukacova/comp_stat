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