# Skrald

#- For data simulation we have created a S3 class called `My_parameters`
#+ The class takes three slots, mu, sigma2 and nu. t. 
#+ The class is created by the function `parameters`
#+ The class has a `simulate` method that generates data from the distribution
#
#
#- For MLE estimation we have created a S3 class called `my_mle`
#+ The class takes two slots, mu_mle and sigma2_mle. 
#+ The class is created by the function `my_mle`
#+ The class has a `plot` that plots the estimated parameters against the true parameters
#
#
#- We test our implementations of estimators using a test function `test_estimator` for $6$ different parameter settings
#+ $(\mu, \sigma^2, \nu) \in \{(0,1,1), (0,1,10), (0,100,15), (10,100,1), (10,3,10), (10,3,15) \}$
#  + For each set of parameters we create an `My_parameters` object and generate data from the distribution. 
#
#```{r}
#set.seed(3564)
#params_list <- list(parameters(mu = 0, sigma2 = 1, nu = 1),
#                    parameters(mu = 0, sigma2 = 1, nu = 10),
#                    parameters(mu = 0, sigma2 = 100, nu = 15),
#                    parameters(mu = 10, sigma2 = 100, nu = 1),
#                    parameters(mu = 10, sigma2 = 3, nu = 10),
#                    parameters(mu = 10, sigma2 = 3, nu = 15))
#```#



#E.logW <- function(x, parameters){
#  mu <- parameters$mu
#  sigma2 <- parameters$sigma2
#  nu <- parameters$nu
#  digamma ((nu + 1) / 2) + log(2) - log(1 + (x - mu)^2 / (nu * sigma2))
#}
#
#Q <- function(X, parameters){
#  mu <- parameters$mu
#  sigma2 <- parameters$sigma2
#  nu <- parameters$nu
#  
#  first_term <- - n * length(X) * log(sqrt(pi * nu * sigma2) * 2^((nu + 1) / 2) * gamma(nu / 2))
#  second_term <- sum(E.logW(X, parameters)) * (nu - 1) / 2
#  third_term <- -sum(E.W(X, parameters) / 2 * (1 + (X - mu)^2 / (nu * sigma2)))
#  
#  first_term + second_term + third_term
#}
