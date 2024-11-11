#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "header_rcpp_functions.h"
#include "sign_func.h"
using namespace Rcpp;
using namespace arma;

// Functions to simulate excess times in situations where the hazard starts out at the out-of-control at a given time point for the linear accelerated time alternative.

double excess_function_to_solve_acc_time_lin_piecewise_out_non_traditional(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    double t,
    double u,
    double linear_predictor,
    double k
) {
  
  double func_eval_at_t = out_of_control_acc_time_lin_cumulative_excess_hazard_add_piecewise_func_single_non_traditional(
    partition_t_vec,
    baseline_vec,
    t,
    k
  );
  
  double final_func_to_solve = 1 - exp(-func_eval_at_t * exp(linear_predictor)) - u; // inverse transform 
  
  return final_func_to_solve;
  
}

double excess_brent_acc_time_lin_piecewise_out_non_traditional(double x1, double x2, double u, double tol, double linear_predictor, double k, arma::vec& partition_t_vec, arma::vec& baseline_vec) { // Implementation of Brent's method from the book "Numerical Recipes: The Art of Scientific Computing"
  
  // Implementation of the Brent algorithm to solve the inverse transform from above. 
  // Adopted from "Numerical Recipes" by Press, Vetterling and Teukolsky (1986).
  
  const int ITMAX = 100;
  // const double EPS = std::numeric_limits<double>::epsilon();
  const double EPS = 3e-8;
  
  double a = x1, b = x2, c = x2, d, e, fa = excess_function_to_solve_acc_time_lin_piecewise_out_non_traditional(partition_t_vec, baseline_vec, a, u, linear_predictor, k), fb = excess_function_to_solve_acc_time_lin_piecewise_out_non_traditional(partition_t_vec, baseline_vec, b, u, linear_predictor, k), fc, p, q, r, s, tol1, xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    throw("Root must be bracketed in brent");
  }
  
  fc = fb;
  
  for (int iter = 0; iter < ITMAX; iter++) {
    
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      
      c = a;
      fc = fa;
      e = b - a;
      d = b - a;
      
    }
    
    if (abs(fc) < abs(fb)) {
      
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;      
      
    }
    
    tol1 = 2.0 * EPS * abs(b) + 0.5 * tol; 
    
    xm = 0.5 * (c - b);
    
    if ((abs(xm) <= tol1 || fb == 0.0)) {
      // std::cout << "Number of iterations: " << iter + 1 << std::endl;
      return b;
      
    }
    
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
      
      s = fb / fa;
      
      if (a == c) {
        
        p = 2.0 * xm * s;
        q = 1.0 - s;
        
      } else {
        
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        
      }
      
      if (p > 0.0) {
        q = -q;
      }
      
      p = abs(p);
      
      double min1 = 3.0 * xm * q - abs(tol1 * q);
      double min2 = abs(e * q);
      
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        
        e = d;
        d = p / q;
        
      } else {
        
        d = xm;
        e = d;
        
      }
      
    } else {
      
      d = xm;
      e = d;
      
    }
    
    a = b;
    fa = fb;
    
    if (abs(d) > tol1) {
      
      b += d;
      
    } else {
      
      b += (abs(tol1) * sign_func(xm));
      
    }
    
    fb = excess_function_to_solve_acc_time_lin_piecewise_out_non_traditional(partition_t_vec, baseline_vec, b, u, linear_predictor, k);
    
  }
  
  throw("Maximum number of iterations exceeded in brent");
  
}

double single_excess_sim_acc_time_lin_piecewise_out_non_traditional(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    double max_follow_up,
    double u,
    double linear_predictor,
    double k
) {
  
  double time_sim = 0;
  
  double func_value_at_end = excess_function_to_solve_acc_time_lin_piecewise_out_non_traditional(
    partition_t_vec, 
    baseline_vec,
    max_follow_up,
    u,
    linear_predictor,
    k
  );
  
  if (func_value_at_end < 0) {
    
    time_sim = max_follow_up; // the inverse transform function to solve never crosses zero --> set the excess time to the max follow up time. 
    
  } else {
    
    time_sim = excess_brent_acc_time_lin_piecewise_out_non_traditional(
      0.0,
      max_follow_up,
      u,
      1e-7,
      linear_predictor, 
      k,
      partition_t_vec, 
      baseline_vec
    );
    
  }
  
  return time_sim;
  
}

//' Function to simulate excess times during a time period where the hazard starts in out-of-control
//' 
//' Function to simulate excess times during a time period where the hazard starts in out-of-control using linear accelerated time alternative
//' 
//' @param partition_t_vec A vector that defines the partition of the time period for the piecewise constant baseline hazard.
//' @param baseline_vec A vector of baseline hazard values defined over the time partition. 
//' @param max_follow_up_vec A vector corresponding to the max follow up times of the observations.  
//' @param u_vec A vector of samples from Unif(0, 1) for inverse transform algorithm. 
//' @param linear_predictor_vec the matrix vector product of the covariate matrix and the regression coefficients. 
//' @param k Linear time acceleration factor for the out-of-control hazard model. 
//' @return time_excess_sim A vector of simulated excess times based on the input of the observations.  
//' @export
// [[Rcpp::export]]
arma::vec vec_excess_sim_acc_time_lin_piecewise_out_non_traditional(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    arma::vec& max_follow_up_vec,
    arma::vec& u_vec,
    arma::vec& linear_predictor_vec,
    double k
) {
  
  int n_iterations = u_vec.size();
  
  arma::vec time_excess_sim(n_iterations);
  
  for (int i = 0; i < n_iterations; i++) {
    
    time_excess_sim(i) = single_excess_sim_acc_time_lin_piecewise_out_non_traditional(
      partition_t_vec,
      baseline_vec,
      max_follow_up_vec(i),
      u_vec(i),
      linear_predictor_vec(i),
      k
    );
    
  }
  
  return time_excess_sim;
  
}