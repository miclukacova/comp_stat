####################################################
########### Piecewise linear Envelope ##############
####################################################

# Slopes
a_i <- function(x_i) df1(x_i) / f_star(x_i) 

# Intercepts
b_i <- function(x_i, a_i) log(f_star(x_i)) - a_i * x_i

# Interval points
z_i <- function(a1, a2, b1, b2) (b2 - b1) / (a1 - a2)

# R_i's
r_i <- function(as, bs, zs, n) {
    1 / as * exp(bs) * (exp(as * zs[2:(n+1)]) - exp(as * zs[1:n]))
}


rf <- function(N, ys) {
  
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
  
  print(R); print(Q)
  
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
    accept[i] <- u[i] <=  f_star(x[i]) / exp(as[int] * x[i] + bs[int])
  }

  return(x[accept])
}

rf(1, c(0.1,0.2,0.3))
