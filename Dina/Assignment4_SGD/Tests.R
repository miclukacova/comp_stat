# SGD tests

#---------------------------- Test of data sampling ----------------------------

N <- 5000
param <- parameters(2,5,1,2)
data <- simulate(param, N)
x <- data$x
y <- data$y


#------------------------------- Test of gradient ------------------------------

gradient(param$par, x, y)





#------------------------------- Test of SGD -----------------------------------

par0 <- c(1,1,1,1)

sgd(par0 = par0, grad = gradient, x = x, y = y, gamma = 0.1, 
    epoch = vanilla, maxiter = 100, sampler = sample)
sgd(par0 = par0, grad = gradient, x = x, y = y, gamma = 0.1, 
    epoch = batch, m = 50, maxiter = 100, sampler = sample)
sgd(par0 = par0, grad = gradient, x = x, y = y, gamma = 0.1, 
    epoch = momentum(), maxiter = 100, sampler = sample)
sgd(par0 = par0, grad = gradient, x = x, y = y, gamma = 0.1,
    m = 100, epoch = adam(), maxiter = 100, sampler = sample, beta1 = 0.5, beta2 = 0.5)



#------------------------------- Test of SGD object-----------------------------

SGD_test <- SGD(par0 = par0, grad = gradient, x = x, y = y, gamma = 0.1, 
                epoch = adam(), maxiter = 100, sampler = sample, cb = SGD_tracer,
                true_par = param$par)

plot(SGD_test, 3)


#--------------------------- Test of vanilla function --------------------------
samp <- sample(N)

vanilla(par = par0, samp = samp, gamma = 0.1, grad = gradient_rcpp, x = x, y = y, N = N)
vanilla_rcpp(par = par0, samp = samp, gamma = 0.1, grad = gradient_rcpp, x = x, y = y, N = N)

microbenchmark::microbenchmark(vanilla(par = par0, samp = samp, gamma = 0.1, grad = gradient_rcpp, x = x, y = y, N = N),
                               vanilla_rcpp(par = par0, samp = samp, gamma = 0.1, grad = gradient_rcpp, x = x, y = y, N = N),
                               times = 10)

#------------------------------- H implementations------------------------------ 

H(x, y, param$par)
H_rcpp(x, y, param$par)

profvis::profvis(sgd(par0 = par0, grad = gradient_rcpp, x = x, y = y, gamma = 0.1, 
    epoch = NULL, maxiter = 100, sampler = sample))



#---------------------------------Test of GD object-----------------------------


GD_test <- GD(x = x, y = y, par0 = c(1,1,1,1), H = H, grad = gradient, d = 0.8, c = 0.1, 
              gamma0 = 1, epsilon = 1e-9, maxiter = 5000, true_par = param$par, clipping = TRUE, momentum = FALSE)
GD_test$est
plot(GD_test, 2)

bench::mark(GD(x = x, y = y, par0 = c(1,1,1,1), H = H, grad = gradient, d = 0.8, c = 0.1, 
              gamma0 = 1, epsilon = 1e-9, maxiter = 5000, true_par = param$par, clipping = TRUE, momentum = FALSE),
            GD(x = x, y = y, par0 = c(1,1,1,1), H = H, grad = gradient, d = 0.8, c = 0.1, 
              gamma0 = 1, epsilon = 1e-9, maxiter = 5000, true_par = param$par, clipping = TRUE, momentum = TRUE),
            times = 10, check = FALSE)
