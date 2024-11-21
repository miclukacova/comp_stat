#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include <vector>
#include <cmath>
#include <numeric>

//------------------------------------SGD algorithm------------------------------------

// [[Rcpp::export]]
NumericVector gradient_rcpp(NumericVector par, NumericVector x, NumericVector y) {
  // Extract parameters
  double alpha0 = par[0];
  double beta0 = par[1];
  double gamma0 = par[2];
  double rho0 = par[3];
  
  // Initialize gradients
  double grad_alpha0 = 0.0;
  double grad_beta0 = 0.0;
  double grad_gamma0 = 0.0;
  double grad_rho0 = 0.0;
  
  int n0 = x.size(); // Number of observations
  
  // Loop over the indices
  for (int i0 = 0; i0 < n0; ++i0) {
    
    // Get individual data point
    double x_i0 = x[i0];
    double y_i0 = y[i0];
    
    // Calculating f(x_i, par)
    double f_x_i0 = gamma0 + (rho0 - gamma0) / (1 + exp(beta0 * log(x_i0) - alpha0));
    
    // Exponential term
    double expbetalogxalpha0 = exp(beta0 * log(x_i0) - alpha0);
    
    // Identical part used in gradients
    double identical_part0 = -2 * (y_i0 - f_x_i0);
    
    // Accumulate gradients for all indices
    grad_alpha0 += (identical_part0 * (rho0 - gamma0) * expbetalogxalpha0) 
      / pow(1 + expbetalogxalpha0, 2);
    grad_beta0 += -(identical_part0 * (rho0 - gamma0) * log(x_i0) * expbetalogxalpha0) 
      / pow(1 + expbetalogxalpha0, 2);
    grad_gamma0 += identical_part0 * (1 - 1 / (1 + expbetalogxalpha0));
    grad_rho0 += identical_part0 / (1 + expbetalogxalpha0);
  }
  
  // Return the mean of accumulated gradients
  return NumericVector::create(grad_alpha0 / n0, grad_beta0 / n0, grad_gamma0 / n0, grad_rho0 / n0);
}


// [[Rcpp::export]]
NumericVector vanilla_rcpp(
    NumericVector par,      // Parameters to update
    IntegerVector samp,     // Indices for sampling
    double gamma,           // Learning rate
    Function grad,          // Gradient function
    NumericVector x,        // x values
    NumericVector y,        // y values
    int N                  // Number of iterations
) {
  
  NumericVector par_new = par;
  
  for (int j = 0; j < N; ++j) {
    int i = samp[j] - 1; // Adjust index for 0-based indexing in C++
    NumericVector grad_i = grad(par, Named("x") = x[i], Named("y") = y[i]);
    
    // Update parameters
    par_new = par - gamma * grad_i;
  }
  
  return par_new;
}


// [[Rcpp::export]]
NumericVector f_rcpp(NumericVector x, NumericVector par) {
  double alpha = par[0];
  double beta = par[1];
  double gamma = par[2];
  double rho = par[3];
  
  int n = x.size();
  NumericVector f(n);
  
  for (int i = 0; i < n; ++i) {
    f[i] = gamma + (rho - gamma) / (1 + std::exp(beta * std::log(x[i]) - alpha));
  }
  
  return f;
}


// [[Rcpp::export]]
double H_rcpp(NumericVector x, NumericVector y, NumericVector par) {
  
  double n = x.size();
  double error_sq = 0.0;
  
  NumericVector f_func = f_rcpp(x, par);
  
  for (int i = 0; i < n; ++i) {
    double error =  y[i] - f_func[i];
    error_sq += error * error;
  }
  
  return error_sq / n;
}



