#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tar_dens_cpp(NumericVector y, NumericVector x, double zx) {
  int n = y.size();
  NumericVector result(n);
  
  // Precompute exp(-sum(exp(y * x)))
  NumericVector exp_yx = exp(y * x);
  double exp_sum = sum(exp_yx);
  
  for (int i = 0; i < n; i++) {
    result[i] = exp(y[i] * zx - exp_sum);
  }
  
  return result;
}

// [[Rcpp::export]]
NumericVector tar_dens_log_difference_cpp(NumericVector y, NumericVector x, double zx) {
  int n = y.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    // Compute the sum of x * exp(y[i] * x)
    double sum_exp = 0;
    for (int j = 0; j < x.size(); j++) {
      sum_exp += x[j] * exp(y[i] * x[j]);
    }
    // Calculate the difference
    result[i] = zx - sum_exp;
  }
  
  return result;
}