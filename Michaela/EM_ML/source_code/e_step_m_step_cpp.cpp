// EM_algorithm.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector E_step_cpp(NumericVector X, NumericVector par, double nu) {
  int n = X.size();
  NumericVector E_W(n);
  double mu = par[0];
  double sigma2 = par[1];
  
  for (int i = 0; i < n; i++) {
    E_W[i] = (nu + 1) / (1 + pow((X[i] - mu), 2) / (nu * sigma2));
  }
  
  return E_W;
}

// [[Rcpp::export]]
NumericVector M_step_cpp(NumericVector E_W, NumericVector X, double nu) {
  int n = X.size();
  double mu_k = 0;
  double sigma2_k = 0;
  double weight_sum = sum(E_W);
  
  // Compute mu_k
  for (int i = 0; i < n; i++) {
    mu_k += E_W[i] * X[i];
  }
  mu_k /= weight_sum;
  
  // Compute sigma2_k
  for (int i = 0; i < n; i++) {
    sigma2_k += E_W[i] * pow((X[i] - mu_k), 2);
  }
  sigma2_k /= (n * nu);
  
  return NumericVector::create(mu_k, sigma2_k);
}

