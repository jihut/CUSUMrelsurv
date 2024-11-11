#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "sign_func.h"
using namespace Rcpp;
using namespace arma;

// Function to simulate excess times based on the smoothed and estimated baseline and cumulative baseline hazard from the semi-parametric model by Perme, Henderson and Stare (2009)
// implemented in relsurv package.
// Corresponds to method = "EM" when using the function rsadd(...).

double rcpp_cumulative_baseline_haz_general_func_EM(Rcpp::Function predict_cumulative_baseline_haz_general_func, double t, double min_cumulative_baseline_event_time) {
  double func_value;
  if (t > min_cumulative_baseline_event_time) {
    NumericVector eval_func = predict_cumulative_baseline_haz_general_func(t);
    func_value = eval_func(0);
  } else { // if the evaluated time point is outside of the time range from the data --> linear interpolation. 
    NumericVector eval_func = predict_cumulative_baseline_haz_general_func(min_cumulative_baseline_event_time);
    func_value = eval_func(0) / min_cumulative_baseline_event_time * t;
  }
  
  return func_value;
  
}

double excess_function_to_solve_EM(Rcpp::Function predict_cumulative_baseline_haz_general_func, double t, double min_cumulative_baseline_event_time, double u, double linear_predictor){
  
  double func_eval_at_t = rcpp_cumulative_baseline_haz_general_func_EM(predict_cumulative_baseline_haz_general_func, t, min_cumulative_baseline_event_time); 
  
  double final_func_to_solve = 1 - exp(-func_eval_at_t * exp(linear_predictor)) - u;  // inverse transform 
  
  return final_func_to_solve; 
  
}

// double sign_func(double num) { // sign function used in the Brent procedure
//   if (num < 0) {
//     return -1;
//   } else if (num > 0) {
//     return 1;
//   } else {
//     return 0;
//   }
// }

double excess_brent_EM(double x1, double x2, double u, double tol, double linear_predictor, double min_cumulative_baseline_event_time, Rcpp::Function predict_cumulative_baseline_haz_general_func) { // Implementation of Brent's method from the book "Numerical Recipes: The Art of Scientific Computing"
  
  // Implementation of the Brent algorithm to solve the inverse transform from above. 
  // Adopted from "Numerical Recipes" by Press, Vetterling and Teukolsky (1986).
  
  const int ITMAX = 100;
  // const double EPS = std::numeric_limits<double>::epsilon();
  const double EPS = 3e-8;
  
  double a = x1, b = x2, c = x2, d, e, fa = excess_function_to_solve_EM(predict_cumulative_baseline_haz_general_func, a, min_cumulative_baseline_event_time, u, linear_predictor), fb = excess_function_to_solve_EM(predict_cumulative_baseline_haz_general_func, b, min_cumulative_baseline_event_time, u, linear_predictor), fc, p, q, r, s, tol1, xm;
  
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
    
    fb = excess_function_to_solve_EM(predict_cumulative_baseline_haz_general_func, b, min_cumulative_baseline_event_time, u, linear_predictor);
    
  }
  
  throw("Maximum number of iterations exceeded in brent");
  
}

double single_excess_sim_EM(
    Rcpp::Function predict_cumulative_baseline_haz_general_func,
    double min_cumulative_baseline_event_time, 
    double u, 
    double linear_predictor,
    double max_follow_up
) {
  
  double time_sim = 0;
  
  double func_value_at_end = excess_function_to_solve_EM(predict_cumulative_baseline_haz_general_func, max_follow_up, min_cumulative_baseline_event_time, u, linear_predictor);
  
  if (func_value_at_end < 0) {
    
    time_sim = max_follow_up; // the inverse transform function to solve never crosses zero --> set the excess time to the max follow up time. 
    
  } else {
    
    time_sim = excess_brent_EM(0.0, max_follow_up, u, 1e-7, linear_predictor, min_cumulative_baseline_event_time, predict_cumulative_baseline_haz_general_func);
    
  }
  
  return time_sim;
  
}

//' Function to simulate excess times based on output from semi-parametric (EM) excess hazard model
//' 
//' Function to simulate excess times based on output from semi-parametric (EM) excess hazard model
//' @param predict_cumulative_baseline_haz_general_func A function that returns predicted value of cumulative baseline based on the smoothed estimated from the EM output from relsurv.
//' @param min_cumulative_baseline_event_time The smallest time to event that was used in the semi-parametric model for estimation.
//' @param u_vec A vector of samples from Unif(0, 1) for inverse transform algorithm. 
//' @param linear_predictor_vec the matrix vector product of the covariate matrix and the regression coefficients. 
//' @param max_follow_up_vec A vector corresponding to the max follow up times of the observations.  
//' @return time_excess_sim A vector of simulated excess times based on the input of the observations.  
//' @export
// [[Rcpp::export]]
arma::vec vec_excess_sim_EM(
    Rcpp::Function predict_cumulative_baseline_haz_general_func,
    double min_cumulative_baseline_event_time,
    arma::vec& u_vec, 
    arma::vec& linear_predictor_vec,
    arma::vec& max_follow_up_vec
) {
  
  int n_iterations = u_vec.size();
  
  arma::vec time_excess_sim(n_iterations);
  
  for (int i = 0; i < n_iterations; i++) {
    
    time_excess_sim(i) = single_excess_sim_EM(predict_cumulative_baseline_haz_general_func, min_cumulative_baseline_event_time, u_vec(i), linear_predictor_vec(i), max_follow_up_vec(i));
    
  }
  
  return time_excess_sim;
  
}