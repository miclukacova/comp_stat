#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List gauss_rej_cpp(int n, double mu, double sigma, double alpha, NumericVector x, double sum1) {
  NumericVector y(n);
  int counter = 0;
  
  for (int i = 0; i < n; ++i) {
    bool reject = true;
    double y0;
    
    while (reject) {
      y0 = R::rnorm(mu, sigma);
      double u = R::runif(0, 1);
      double density = R::dnorm(y0, mu, sigma, false);

      // Compute f^*(y0) directly in C++
      double term1 = y0 * sum1;
      double sum2 = 0.0;
      for (int j = 0; j < x.size(); ++j) {
        sum2 += exp(y0 * x[j]);
      }
      
      double f_value = exp(term1 - sum2);
      reject = u > f_value / density * alpha;
      counter++;
    }
    y[i] = y0;
  }
  
  return List::create(Named("sample") = y,
                      Named("accept") = static_cast<double>(n) / counter);
}

