#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradient_rcpp(NumericVector par, NumericVector x, NumericVector y) {
  double alpha0 = par[0];
  double beta0 = par[1];
  double gamma0 = par[2];
  double rho0 = par[3];
  
  double grad_alpha0 = 0.0;
  double grad_beta0 = 0.0;
  double grad_gamma0 = 0.0;
  double grad_rho0 = 0.0;
  
  int n0 = x.size();
  for (int i0 = 0; i0 < n0; ++i0) {
    double x_i0 = x[i0];
    double y_i0 = y[i0];
    double f_x_i0 = gamma0 + (rho0 - gamma0) / (1 + exp(beta0 * log(x_i0) - alpha0));
    double expbetalogxalpha0 = exp(beta0 * log(x_i0) - alpha0);
    double identical_part0 = -2 * (y_i0 - f_x_i0);
    
    grad_alpha0 += (identical_part0 * (rho0 - gamma0) * expbetalogxalpha0) 
      / pow(1 + expbetalogxalpha0, 2);
    grad_beta0 += -(identical_part0 * (rho0 - gamma0) * log(x_i0) * expbetalogxalpha0) 
      / pow(1 + expbetalogxalpha0, 2);
    grad_gamma0 += identical_part0 * (1 - 1 / (1 + expbetalogxalpha0));
    grad_rho0 += identical_part0 / (1 + expbetalogxalpha0);
  }
  
  return NumericVector::create(grad_alpha0 / n0, grad_beta0 / n0, grad_gamma0 / n0, grad_rho0 / n0);
}

// [[Rcpp::export]]
NumericVector vanilla_rcpp(
    NumericVector par,
    IntegerVector samp,
    double gamma,
    Function grad,
    NumericVector x,
    NumericVector y,
    int N
) {
  NumericVector par_new = clone(par);
  for (int j = 0; j < N; ++j) {
    int i = samp[j] - 1;
    NumericVector grad_i = grad(par_new, Named("x") = x[i], Named("y") = y[i]);
    par_new = par_new - gamma * grad_i;
  }
  return par_new;
}

// [[Rcpp::export]]
List SGD_cpp(
    NumericVector par0,
    NumericVector x,
    NumericVector y,
    Function grad,
    NumericVector gamma,
    int maxiter = 100,
    double epsilon = 1e-6
) {
  int N = x.size();
  NumericVector par = clone(par0);
  NumericVector par_prev = clone(par0);
  
  for (int iter = 0; iter < maxiter; ++iter) {
    IntegerVector samp = sample(N, N, false);  // Random sampling for SGD
    
    par = vanilla_rcpp(par, samp, gamma[iter % gamma.size()], grad, x, y, N);
    
    double diff_norm = sum(pow(par - par_prev, 2));
    if (diff_norm <= epsilon * (sum(pow(par, 2)) + epsilon)) {
      break;
    }
    par_prev = clone(par);
  }
  
  return List::create(
    Named("est") = par,
    Named("iterations") = maxiter
  );
}
