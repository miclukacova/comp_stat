
#------------------------------------------------------------------------------#
#           Poisson target distribution and log and derivative hereof          #
#------------------------------------------------------------------------------#

# Target distribution
target_distribution_pois <- function(y){
  x <- poisson$x
  z <- poisson$z
  return(prod(exp(y*z*x - exp(y*x))))
}
target_distribution_pois <- Vectorize(target_distribution_pois, vectorize.args = "y")

# Log target distribution
log_target_distribution_pois <- function(y){
  x <- poisson$x
  z <- poisson$z
  return(sum(y*z*x - exp(y*x)))
}
log_target_distribution_pois <- Vectorize(log_target_distribution_pois, vectorize.args = "y")

# Derivative of log target distribution
log_target_distribution_diff_pois <- function(y){
  x <- poisson$x
  z <- poisson$z
  return(sum(z*x - x*exp(y*x)))
}
log_target_distribution_diff_pois <- Vectorize(log_target_distribution_diff_pois, vectorize.args = "y")







target_distribution_pois1 <- function(y, x, z, xz){
  return(exp(sum(y*xz - exp(y*x))))
}
target_distribution_pois1 <- Vectorize(target_distribution_pois1, vectorize.args = "y")


target_distribution_pois2 <- function(y, x, z, xz){
  return(exp(rowSums(y %*% t(xz) - exp(y %*% t(x)))))
}

target_distribution_pois3 <- function(y, x, z, xz){
  return(exp(colSums(xz %*% t(y) - exp(x %*% t(y)))))
}




























#------------------------------------------------------------------------------#
#                   Functions to help with different calculations              #
#------------------------------------------------------------------------------#

# Proportion between proposal and target distribution
proportion_function <- function(target_distribution, proposal_distribution, y, ...){
  return(proposal_distribution(y, ...)/target_distribution(y))
}

log_proportion_function_pois <- function(y, sigma, mu){
  d <- log(sqrt(2*pi)*sigma)
  proposal <- d * 0.5*(y - mu)^2/sigma^2
  target <- sum(y*poisson$z*poisson$x - exp(y*poisson$x))
  
  return(proposal - target)
}
log_proportion_function_pois <- Vectorize(log_proportion_function_pois, vectorize.args = "y")



alpha_star_sigma <- function(sigma, y_vec){
  
  # prop <- log_proportion_function_pois(y = y_vec, sigma = sigma, mu = mu_hat)
  
  prop <- dnorm(y_vec, mean = mu_hat, sd = sigma) / target_distribution_pois(y_vec)
  
  alpha_star <- min(prop[prop > 0])
  
  return(alpha_star)
}

sigma_hat_func <- function(sigma_seq, y_vec){
  
  alpha_stars <- sapply(sigma_seq, FUN = alpha_star_sigma, y_vec = y_vec)
  
  index <- which.max(alpha_stars)
  
  sigma_hat <- sigma_seq[index]
  
  return(sigma_hat)
}


# alpha_prime_sigma <- function(sigma, x_seq) {
#   
#   alpha_prime <- optimize(f = function(y) 
#     log_proportion_function_pois(y = y,
#                         mu = mu_hat, 
#                         sigma = sigma),
#     interval = c(0, 1), maximum = FALSE)$objective
#   
#   return(alpha_prime)
# }




# z function
z_function <- function(a,b,min,max){
  
  z <- (b[-1] - b[-length(b)])/(a[-length(a)] - a[-1])
  
  z <- c(min, z, max)
  return(z)
}


# H_i function for adaptive envelope
H_i <- function(a, b, z, break_points){
  
  z_0 <- z[-length(z)]
  z_1 <- z[-1]
  
  return(1 / a * exp(b) * (exp(a*z_1) - exp(a*z_0)))
  
}

H_i_inv <- function(a, b, z, Q_is, dq){
  
}




dq_in_Q_i <- function(dq, Q_is, z){
  # Note that the addition of min(z) and is slightly iffy, so this might need to be changed if things look weird
  lower <- outer(c(0, Q_is[-length(Q_is)]), dq, "<") # Note that the order here is important for the dimension of the matrix
  upper <- outer(Q_is, dq, ">=")
  
  return(lower * upper)
}




