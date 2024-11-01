#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec EM_alg_arma(const arma::vec& x, Rcpp::List param, double ny, int max_iter = 20, double epsilon = 1e-10) {
  // Initialize parameters
  double mu_mark = Rcpp::as<double>(param["mu"]);
  double sigma_mark = Rcpp::as<double>(param["sigma"]);
  double k = (ny + 1) / 2.0;  // Define k based on ny
  int n = x.n_elem;           // Size of the vector x
  int i = 0;                  // Iteration counter
  
  while (i < max_iter) {
    // Store old parameters
    double mu_old = mu_mark;
    double sigma_old = sigma_mark;
    
    // E-step: Calculate t_old (using element-wise operations)
    arma::vec t_old = 2 / (1 + arma::pow(x - mu_old, 2) / (ny * sigma_old));
    
    // M-step: Update mu and sigma
    mu_mark = arma::dot(t_old, x) / arma::sum(t_old);  // Weighted sum for mu
    sigma_mark = (1.0 / (n * ny)) * k * arma::dot(t_old, arma::pow(x - mu_mark, 2));
    
    // Convergence check
    if (std::pow(mu_mark - mu_old, 2) + std::pow(sigma_mark - sigma_old, 2) <= 
        epsilon * (std::pow(mu_mark, 2) + std::pow(sigma_mark, 2) + epsilon)) {
      break;
    }
    
    i++;
  }
  
  // Return final values for mu and sigma
  return arma::vec({mu_mark, sigma_mark});
}