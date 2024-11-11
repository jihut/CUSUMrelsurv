#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include "sign_func.h"
using namespace Rcpp;
using namespace arma;

// Function to simulate excess times from a piecewise baseline excess hazard.

int find_interval_single(double& x, arma::vec& vec) {
  int final_index = 0;
  
  arma::uvec indices = arma::find(vec <= x);
  if (!indices.is_empty()) {
    final_index = indices.max();
  }
  
  return final_index; 
  // flaw if this is used for general purposes as it will return 0 for observations less than the first element too, but this is needed to define uvec for later applications.
  // the first value of vec will always be 0 in our case and there will not be any negative follow up time.
}

double cumulative_baseline_excess_hazard_piecewise_func_single(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    double& t
) {
  vec interval_length_partition_t_vec = diff(partition_t_vec); // time between two time points from the input time vector
  int num_interval_length_partition_t_vec = interval_length_partition_t_vec.size(); // number of intervals between the time points given
  uvec index_vec = linspace<uvec>(0, num_interval_length_partition_t_vec - 1, num_interval_length_partition_t_vec);
  vec new_haz_vec = baseline_vec.elem(index_vec); // relevant hazard values during these intervals
  vec pre_cumulative_baseline_vec = cumsum(new_haz_vec % interval_length_partition_t_vec); // compute cumulative hazard at these time points
  vec front_cumulative_baseline_vec = {0};
  vec cumulative_baseline_vec = join_cols(front_cumulative_baseline_vec, pre_cumulative_baseline_vec);
  int index = find_interval_single(t, partition_t_vec);
  double relevant_cumulative_baseline_excess_hazard_piecewise = cumulative_baseline_vec(index) + (t - partition_t_vec(index)) * baseline_vec(index);
  
  return relevant_cumulative_baseline_excess_hazard_piecewise;
}

double excess_function_to_solve_piecewise(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    double t,
    double u,
    double linear_predictor
) {
  
  double func_eval_at_t = cumulative_baseline_excess_hazard_piecewise_func_single(
    partition_t_vec,
    baseline_vec,
    t
  );
  
  double final_func_to_solve = 1 - exp(-func_eval_at_t * exp(linear_predictor)) - u; // inverse transform
  
  return final_func_to_solve;
  
}

double excess_brent_piecewise(double x1, double x2, double u, double tol, double linear_predictor, arma::vec& partition_t_vec, arma::vec& baseline_vec) { // Implementation of Brent's method from the book "Numerical Recipes: The Art of Scientific Computing"
  
  // Implementation of the Brent algorithm to solve the inverse transform from above. 
  // Adopted from "Numerical Recipes" by Press, Vetterling and Teukolsky (1986).
  
  const int ITMAX = 100;
  // const double EPS = std::numeric_limits<double>::epsilon();
  const double EPS = 3e-8;
  
  double a = x1, b = x2, c = x2, d, e, fa = excess_function_to_solve_piecewise(partition_t_vec, baseline_vec, a, u, linear_predictor), fb = excess_function_to_solve_piecewise(partition_t_vec, baseline_vec, b, u, linear_predictor), fc, p, q, r, s, tol1, xm;
  
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
    
    fb = excess_function_to_solve_piecewise(partition_t_vec, baseline_vec, b, u, linear_predictor);
    
  }
  
  throw("Maximum number of iterations exceeded in brent");
  
}

double single_excess_sim_piecewise(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    double max_follow_up,
    double u,
    double linear_predictor
) {
  
  double time_sim = 0;
  
  double func_value_at_end = excess_function_to_solve_piecewise(
    partition_t_vec, 
    baseline_vec,
    max_follow_up,
    u,
    linear_predictor
  );
  
  if (func_value_at_end < 0) {
    
    time_sim = max_follow_up; // the inverse transform function to solve never crosses zero --> set the excess time to the max follow up time. 
    
  } else {
    
    time_sim = excess_brent_piecewise(
      0.0,
      max_follow_up,
      u,
      1e-7,
      linear_predictor, 
      partition_t_vec, 
      baseline_vec
    );
    
  }
  
  return time_sim;
  
}

//' Function to simulate excess times based on output from piecewise constant baseline excess hazard model
//' 
//' Function to simulate excess times based on output from piecewise constant baseline excess hazard model
//' @param partition_t_vec A vector that defines the partition of the time period for the piecewise constant baseline hazard.
//' @param baseline_vec A vector of baseline hazard values defined over the time partition. 
//' @param u_vec A vector of samples from Unif(0, 1) for inverse transform algorithm. 
//' @param linear_predictor_vec the matrix vector product of the covariate matrix and the regression coefficients. 
//' @param max_follow_up_vec A vector corresponding to the max follow up times of the observations.  
//' @return time_excess_sim A vector of simulated excess times based on the input of the observations.  
//' @export
// [[Rcpp::export]]
arma::vec vec_excess_sim_piecewise(
    arma::vec& partition_t_vec, 
    arma::vec& baseline_vec,
    arma::vec& max_follow_up_vec,
    arma::vec& u_vec,
    arma::vec& linear_predictor_vec
) {
  
  int n_iterations = u_vec.size();
  
  arma::vec time_excess_sim(n_iterations);
  
  for (int i = 0; i < n_iterations; i++) {
    
    time_excess_sim(i) = single_excess_sim_piecewise(
      partition_t_vec,
      baseline_vec,
      max_follow_up_vec(i),
      u_vec(i),
      linear_predictor_vec(i)
    );
    
  }
  
  return time_excess_sim;
  
}