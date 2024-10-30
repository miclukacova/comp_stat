library(parallel)

piece_lin_rejec_samp_par <- function(N, ys) {
  
  # Calculating a's, b's, z's
  as <- sapply(ys, a_i, simplify = TRUE)
  bs <- mapply(FUN = b_i, ys, as)
  n <- length(bs)
  zs <- c(-Inf, mapply(FUN = z_i, as[1:(n-1)], as[2:n], bs[1:(n-1)], bs[2:n]), Inf)
  
  # Bookkeeping: I_i integrals
  R <- r_i(as, bs, zs, n)
  Q <- c(0, cumsum(R))
  
  # Drawing from piecewise linear density and uniform
  u0 <- Q[n + 1] * runif(N)
  u <- runif(N)
  
  # Determine the interval that each point belongs to
  geq_z <- outer(u0, Q[1:n], FUN = function(y1, y2) y1 > y2)
  leq_z <- outer(u0, Q[2:(n+1)], FUN = function(y1, y2) y1 <= y2)
  I <- geq_z & leq_z
  
  # Parallelized sampling
  results <- mclapply(1:N, function(i) {
    int <- which(I[i,] == 1)
    x_i <- log((u0[i] - Q[int]) * as[int] * exp(-bs[int]) + exp(as[int] * zs[int])) / as[int]
    accept <- u[i] <= f_star(x_i) / exp(as[int] * x_i + bs[int])
    if (accept) x_i else NA  # Return x_i if accepted, else NA
  }, mc.cores = detectCores()) # Adjust cores if needed
  
  # Filter accepted values and remove NAs
  x <- unlist(results)
  x <- x[!is.na(x)]
  
  return(x)
}
