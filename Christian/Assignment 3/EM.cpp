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


//------------------------------------EM algorithm------------------------------------

// [[Rcpp::export]]
double alpha_prime_func_cpp(double nu){
  return (nu + 1) / 2.0;
}

// [[Rcpp::export]]
NumericVector beta_prime_func_cpp(NumericVector x, double nu, double mu, double sigma) {
  int n = x.size();
  NumericVector beta_prime(n);
  
  for (int i = 0; i < n; i++) {
    beta_prime[i] = 0.5 * (1 + std::pow(x[i] - mu, 2) / (nu * std::pow(sigma, 2)));
  }
  
  return beta_prime;
}


// [[Rcpp::export]]
double mu_hat_func_cpp(NumericVector x, NumericVector beta_prime) {
  int n = x.size();
  double sum_x_beta = 0.0;
  double sum_beta_inv = 0.0;
  
  for (int i = 0; i < n; i++) {
    sum_x_beta += x[i] / beta_prime[i];
    sum_beta_inv += 1.0 / beta_prime[i];
  }
  
  // mu_hat = sum(x / beta_prime) / sum(1 / beta_prime)
  return sum_x_beta / sum_beta_inv;
}

// [[Rcpp::export]]
double sigma_hat_func_cpp(NumericVector x, double nu, double mu_hat, 
                          double alpha_prime, NumericVector beta_prime) {
  int n = x.size();
  double sum_squared_diff = 0.0;
  
  for (int i = 0; i < n; i++) {
    sum_squared_diff += std::pow(x[i] - mu_hat, 2) / beta_prime[i];
  }
  
  // sigma_hat = sqrt(alpha_prime / nu * mean((x - mu_hat)^2 / beta_prime))
  return std::sqrt(alpha_prime / nu * (sum_squared_diff / n));
}


// [[Rcpp::export]]
NumericVector EM_algorithm_full_cpp(NumericVector x, NumericVector par0, double nu, 
                               int maxiter = 1000, double tolerance = 1e-6) {
  
  double mu_hat = par0[0];
  double sigma_hat = par0[1];
  double alpha_prime = alpha_prime_func_cpp(nu);
  
  
  for (int i = 0; i < maxiter; i++) {
    double mu_old = mu_hat;
    double sigma_old = sigma_hat;
    
    // E-step: Calculate beta_prime based on current mu and sigma estimates
    NumericVector beta_prime = beta_prime_func_cpp(x, nu, mu_old, sigma_old);

    // M-step: Update mu_hat and sigma_hat
    mu_hat = mu_hat_func_cpp(x, beta_prime);
    sigma_hat = sigma_hat_func_cpp(x, nu, mu_hat, alpha_prime, beta_prime);
    
    // Check convergence criteria
    if (std::abs(mu_hat - mu_old) < tolerance && std::abs(sigma_hat - sigma_old) < tolerance) {
      break;
    }
  }
  
  return NumericVector::create(mu_hat, sigma_hat);
}









//---------------------------------Gradient descent---------------------------------

// [[Rcpp::export]]

double marginal_log_likelihood_cpp(NumericVector x, double mu, double sigma, double nu){
  int n = x.size();
  double sum_log = 0.0;
  
  for (int i = 0; i < n; i++) {
    sum_log += std::log(1 + std::pow(x[i] - mu, 2) / (nu * std::pow(sigma, 2)));
  }
  
  return -n * std::log(sigma) - (nu + 1) / 2 * sum_log;
}





// [[Rcpp::export]]
NumericVector d_marginal_log_likelihood_cpp(NumericVector x, double mu, double sigma, double nu){
  int n = x.size();
  double d_sigma = -n / sigma;
  double d_mu = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double diff = x[i] - mu;
    double diff_sq = diff * diff;
    double denom1 = sigma * sigma * sigma * (1 + diff_sq / (nu * sigma * sigma));
    double denom2 = sigma * sigma + diff_sq / nu;
    d_sigma += (nu + 1) / nu * diff_sq / denom1;
    d_mu += (nu + 1) / nu * diff / denom2;
  }
  
  return NumericVector::create(d_mu, d_sigma);
}



// [[Rcpp::export]]
NumericVector grad_clipped_cpp(NumericVector grad, double grad_norm, double threshold = 1.0) {
  if (grad_norm > threshold) {
    grad = grad / grad_norm * threshold;
  }
  return grad;
}




