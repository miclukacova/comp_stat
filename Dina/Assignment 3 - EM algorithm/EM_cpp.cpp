#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector EM_alg_cpp(NumericVector x, List param, double ny, int max_iter = 20, double epsilon = 1e-10) {
  // Get initial parameters
  double mu_mark = as<double>(param["mu"]);
  double sigma_mark = as<double>(param["sigma"]);
  double k = (ny + 1) / 2.0; // Calculating k based on ny
  int n = x.size(); // Length of vector x
  int i = 0; // Iteration counter
  
  while (i < max_iter) {
    // Save old parameters
    double mu_old = mu_mark;
    double sigma_old = sigma_mark;
    
    // E-step: Calculate t_old
    NumericVector t_old = 2 / (1 + pow(x - mu_old, 2) / (ny * sigma_old));
    
    // M-step: Update mu and sigma
    mu_mark = sum(t_old * x) / sum(t_old);
    sigma_mark = (1.0 / (n * ny)) * k * sum(t_old * pow(x - mu_mark, 2));
    
    // Check convergence
    if (pow(mu_mark - mu_old, 2) + pow(sigma_mark - sigma_old, 2) <= 
        epsilon * (pow(mu_mark, 2) + pow(sigma_mark, 2) + epsilon)) {
      break;
    }
    
    i++;
  }
  
  return NumericVector::create(mu_mark, sigma_mark);
}
