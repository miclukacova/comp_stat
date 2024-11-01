### The profiling slide:
---
  ###CV Bandwidth selection
  By construction, the CV bandwidth selection is more computationally expensive. So we investigate whether it is possible to optimize it. First we profile to identify the bottlenecks
```{r, echo = FALSE}
profvis({
  cv_bw_M <- function(x, k = 10, h = seq(0.03, 10, by = 0.01)){
    #Partition
    n <-length(x)
    groups <- sample(rep(1:k, length.out = n))
    
    #Number of observations in each group:
    N <- sapply(seq(1,k), function(x) sum (groups == x), simplify = T)
    
    #bandwidth log likelihoods:
    bw_l <- numeric(length(h))
    
    for (j in seq_along(h)) {
      #kernel density estimates:
      f_hat_i <- numeric(length(x))
      
      #Calculating CV density estimates
      for (i in seq_along(x)) {
        k <- groups[i]
        x_m_i <- x[groups != k]
        condition <- abs((x[i] - x_m_i)/h[j]) <= 1
        f_hat_i[i] <- 3 / (N[k] * 4 * h[j]) * sum((1-(x[i] - x_m_i)^2/h[j]^2) * condition)
        #if(any(f_hat_i < 0 )) {print(f_hat_i)} #Check for negative values of f_hat
      }
      
      #Calculating log likelihood
      bw_l[j] <- sum(log(f_hat_i))
      #print(c(bw_l[j], h[j]))
      
    }
  }
  cv_bw_M(x = rnorm(1000, 0, 1), k = 10, h = seq(0.03, 10, by = 0.01))
})
```
