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