#include <Rcpp.h>
using namespace Rcpp;

// Example code fragment translated to Rcpp
// Assumes kern_vals_h is a precomputed matrix, and N is a vector of counts per group

// [[Rcpp::export]]
NumericVector inner_loop_rcpp(NumericVector x, double h, IntegerVector groups, 
                                NumericMatrix kern_vals_h, IntegerVector N) {
  int n = x.size();
  NumericVector f_hat_i(n);
  
  for (int i = 0; i < n; ++i) {
    int kk = groups[i]; // Group for cross-validation
    double sum_vals = 0.0;
    
    // Compute the condition for each element
    for (int m = 0; m < n; ++m) {
      if (m != i && std::abs(x[m] - x[i]) <= h && groups[m] != kk) {
        sum_vals += kern_vals_h(i, m);
      }
    }
    
    // Calculate density estimate for each observation
    f_hat_i[i] = (N[kk - 1] > 0) ? (sum_vals / N[kk - 1]) : 0;
  }
  
  f_hat_i = 3 * f_hat_i / (4 * h);
  
  return f_hat_i;
}

