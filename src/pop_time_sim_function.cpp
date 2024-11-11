#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <string.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace Rcpp;
using namespace arma;
#include "header_rcpp_functions.h"

// Functions to simulate population times from population tables.

NumericVector rcpp_seq(double from_input, double to_input, double by_input = 1.0) { // faster with NumericVector here, no problem when using vec later without converting too.
  
  int n = (to_input - from_input) / by_input + 1; // check how many elements the sequence will have
  
  NumericVector seq(n); // initialize a vector with length gievn above
  
  for (int i = 0; i < n; i++) { // loop to get the sequence of interest
    seq[i] = from_input + i * by_input;
  }
  
  return seq;
  
}

double submat(arma::mat& X, int age, int year) { // no need to export this later on
  
  double haz_value = 0.0;
  int n = X.n_rows;
  
  for (int i = 0; i < n; i++) {
    if(X(i, 0) == year && X(i, 1) == age) { // scan over rows in X to find the row corresponding to the specified combination of year and age
      haz_value = X(i, 2); // second column corresponds to the hazard value
      break;
    }
  }
  
  return haz_value;
}

arma::mat haz_output_single(double age, double start_year, double end_year, double arrival_time, double table_end_year, arma::mat& pop_matrix){
  
  double step_age = ceil(age) - age; // how long until the hazard changes due to age change
  double step_year = ceil(start_year + arrival_time) - (start_year + arrival_time); // how long until the hazard changes due to year change
  double max_follow_up = end_year - arrival_time - start_year; // max follow up time
  
  NumericVector final_step_age_vec = rcpp_seq(step_age, max_follow_up, 1); // all time points (where arrival equals to t = 0) where hazard change due to age occurs until max follow up
  NumericVector final_step_year_vec = rcpp_seq(step_year, max_follow_up, 1);
  
  int length_of_step_age = final_step_age_vec.size();
  int length_of_step_year = final_step_year_vec.size();
  
  arma::vec time_vec;
  
  if (step_age == 0 || step_year == 0) {
    
    arma::vec test_time_vec(length_of_step_age + length_of_step_year); // vector of time points where hazard changes - length_of_step_age + length_of_step_year elements since t = 0 already is included
    
    for (int i = 0; i < length_of_step_age; ++i){
      test_time_vec[i] = final_step_age_vec[i];
    }
    
    for(int i = length_of_step_age; i < (length_of_step_age + length_of_step_year); ++i){
      test_time_vec[i] = final_step_year_vec[i - length_of_step_age];
    }
    
    time_vec = test_time_vec;
    
  } else {
    
    arma::vec test_time_vec(length_of_step_age + length_of_step_year + 1); // length_of_step_age + length_of_step_year + 1 since we need the first element to be at t = 0
    
    for (int i = 1; i < (length_of_step_age + 1); ++i) {
      test_time_vec[i] = final_step_age_vec[i - 1];
    }
    
    for (int i = length_of_step_age + 1; i < (length_of_step_age + length_of_step_year + 1); ++i) {
      test_time_vec[i] = final_step_year_vec[i - (length_of_step_age + 1)];
    }
    
    time_vec = test_time_vec;
    
  }
  arma::uvec index_time_before_max_follow_up = find(time_vec < max_follow_up);
  arma::vec new_time_vec = time_vec.elem(index_time_before_max_follow_up); // only need the time points before max_follow_up
  if (step_age == step_year) {
    new_time_vec = unique(new_time_vec); // duplicate when step_age == step_year
  }
  std::sort(new_time_vec.begin(), new_time_vec.end());
  
  arma::vec age_vec_at_different_times = floor(age + new_time_vec); // age during the time points given above
  arma::vec year_vec_at_different_times = floor(start_year + arrival_time + new_time_vec); // year during the time points given above
  arma::uvec limit_age_indices = find(age_vec_at_different_times > 110);
  age_vec_at_different_times.elem(limit_age_indices).fill(110); // max value at 110 based on the formatted population table
  arma::uvec limit_table_year_indices = find(year_vec_at_different_times >= table_end_year + 1);
  year_vec_at_different_times.elem(limit_table_year_indices).fill(table_end_year); // in case we have a year value outside the range in population table
  
  arma::vec haz_values(new_time_vec.size()); // storage for the hazard values based on the combination of year and age

  for (int j = 0; j < new_time_vec.size(); ++j) {
    haz_values[j] = submat(pop_matrix, age_vec_at_different_times[j], year_vec_at_different_times[j]); // the last two arguments of submat_rcpp should be integers, will automatically round down.
  }
  
  arma::mat result = join_rows(new_time_vec, haz_values);
  
  return result;
}

arma::vec pop_haz_func(arma::vec& age_vec, arma::vec& gender_vec, arma::vec& arrival_time_vec, arma::vec& t_vec, arma::mat& pop_matrix_male, arma::mat& pop_matrix_female, double start_year, double end_year, double table_end_year) {
  
  int n_elements = age_vec.size();
  
  arma::vec result(n_elements);
  
  for (int i = 0; i < n_elements; i++) {
    
    arma::mat pop_matrix; 
    
    if (gender_vec[i] == 1) {
      pop_matrix = pop_matrix_male;
    } else if (gender_vec[i] == 2) { 
      pop_matrix = pop_matrix_female;
    }
    
    arma::mat haz_matrix = haz_output_single(age_vec[i], start_year, end_year, arrival_time_vec[i], table_end_year, pop_matrix); // uses the function above to choose out the correct row from population table for the given observation.
    
    arma::vec haz_time_intervals = haz_matrix.col(0);
    arma::vec haz_values = haz_matrix.col(1);
    
    int relevant_index = find_interval_single(t_vec(i), haz_time_intervals);
    
    result(i) = haz_values(relevant_index);
    
  }
  
  return result;
  
}

double cumulative_pop_function(double t, arma::vec& time_vec, arma::vec& haz_vec){ // input like time_vec and haz_vec should be the output from haz_output_single
  
  double cumulative_value;
  
  if (time_vec.size() != 1) {
    
    // Want to construct a vector that gives the cumulative hazard at the time points where the hazard changes
    
    arma::vec interval_length_partition_t_vec = diff(time_vec); // time between two time points from the input time vector
    int num_interval_length_partition_t_vec = interval_length_partition_t_vec.size(); // number of intervals between the time points given
    arma::uvec index_vec = linspace<uvec>(0, num_interval_length_partition_t_vec - 1, num_interval_length_partition_t_vec);
    arma::vec new_haz_vec = haz_vec.elem(index_vec); // relevant hazard values during these intervals
    
    arma::vec pre_cumulative_baseline_vec = cumsum(new_haz_vec % interval_length_partition_t_vec); // compute cumulative hazard at these time points
    arma::vec front_cumulative_baseline_vec = {0};
    arma::vec cumulative_baseline_vec = join_cols(front_cumulative_baseline_vec, pre_cumulative_baseline_vec);
    
    int index = 0;
    
    if (t < time_vec[time_vec.size() - 1]) { // find which time interval between two time points t belongs in if t is less than the largest time point given
      
      for (int i = 0; i < time_vec.size(); i++) {
        
        if (time_vec[i] <= t && t < time_vec[i + 1]) {
          break;
        } else {
          index += 1;
        }
        
      }
      
    } else {
      
      index = time_vec.size() - 1;
      
    }
    
    cumulative_value = cumulative_baseline_vec[index] + (t - time_vec[index]) * haz_vec[index]; // final cumulative hazard value making use of the computation at the time points given as input
    
  } else {
    
    cumulative_value = haz_vec[0] * t;
    
  }
  
  return cumulative_value;
  
}

double pop_function_to_solve(double t, double u, arma::vec& time_vec, arma::vec& haz_vec){
  return 1 - exp(-cumulative_pop_function(t, time_vec, haz_vec)) - u; // comes from the inverse transform method, need to solve/find the root of this equation with respect to t to get the simulated population time
}

double sign_func(double num) { // sign function used in the Brent procedure
  if (num < 0) {
    return -1;
  } else if (num > 0) {
    return 1;
  } else {
    return 0;
  }
}

double pop_brent(double x1, double x2, double u, double tol, vec& time_vec, vec& haz_vec) { // Implementation of Brent's method from the book "Numerical Recipes: The Art of Scientific Computing"
  
  // Implementation of the Brent algorithm to solve the inverse transform from above. 
  // Adopted from "Numerical Recipes" by Press, Vetterling and Teukolsky (1986).
  
  const int ITMAX = 100;
  // const double EPS = std::numeric_limits<double>::epsilon();
  const double EPS = 3e-8;
  
  double a = x1, b = x2, c = x2, d, e, fa = pop_function_to_solve(a, u, time_vec, haz_vec), fb = pop_function_to_solve(b, u, time_vec, haz_vec), fc, p, q, r, s, tol1, xm;
  
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
    
    fb = pop_function_to_solve(b, u, time_vec, haz_vec);
    
  }
  
  throw("Maximum number of iterations exceeded in brent");
  
}

arma::vec single_pop_sim(
    double age,
    int gender,
    double start_year,
    double end_year,
    double table_end_year,
    double arrival_time,
    arma::mat& pop_matrix_male,
    arma::mat& pop_matrix_female,
    double u
) {
  
  double max_follow_up_time = end_year - arrival_time - start_year;
  
  arma::mat haz_matrix;
  if (gender == 1) {
    haz_matrix = haz_output_single(age, start_year, end_year, arrival_time, table_end_year, pop_matrix_male);
  } else if (gender == 2){
    haz_matrix = haz_output_single(age, start_year, end_year, arrival_time, table_end_year, pop_matrix_female);
  }
  arma::vec time_vec = haz_matrix.col(0);
  arma::vec haz_vec = haz_matrix.col(1);
  
  double time_sim = 0;
  double delta_p = 0; // initialize event indicator
  double func_value_end = pop_function_to_solve(max_follow_up_time, u, time_vec, haz_vec); // the value of the function we want to find the root of at the maximum follow up the individual can have
  if (func_value_end < 0) { // if the value is less than 0 --> the function will not cross zero during the possible follow up interval --> set the simulated population time to the maximum follow up time
    time_sim = max_follow_up_time;
  } else {
    time_sim = pop_brent(0.0, max_follow_up_time, u, 1e-7, time_vec, haz_vec);
    delta_p = 1;
  }
  
  return {time_sim, delta_p};
}

//' Function to simulate population times based on life tables
//' 
//' Function to simulate population times based on life tables
//' @param age_vec A vector containing age of the observations when they arrive in the study. 
//' @param gender_vec A vector containing gender indicator - 0 for male, 1 for female. 
//' @param start_year Calendar year of when the monitoring starts. 
//' @param end_year Calendar year of when the monitoring ends. 
//' @param table_end_year The last calendar year appearing in the population table. 
//' @param arrival_time_vec A vector containing the amount of time units after start year the observations arrive.
//' @param pop_matrix_male Population table for male, should be a matrix with first column corresponding to year, second column to age and third column to the population hazard value. 
//' @param pop_matrix_female Population table for female, should be a matrix with first column corresponding to year, second column to age and third column to the population hazard value. 
//' @param u A vector of samples from Unif(0, 1) for inverse transform algorithm. 
//' @return A matrix where first column corresponds to the population time and the second is the event indicator. 
//' @export
// [[Rcpp::export]]
arma::mat vec_pop_sim(arma::vec& age_vec,
                     arma::vec& gender_vec,
                     double& start_year,
                     double& end_year,
                     double& table_end_year,
                     arma::vec& arrival_time_vec,
                     arma::mat& pop_matrix_male,
                     arma::mat& pop_matrix_female,
                     arma::vec& u){ // final function to simulate the population times
  
  int n_iterations = age_vec.size();
  arma::vec time_pop_sim(n_iterations);
  arma::vec delta_pop_sim(n_iterations);
  
  for (int i = 0; i < n_iterations; i++) {
    
    vec pop_sim = single_pop_sim(age_vec(i), gender_vec(i), start_year, end_year, table_end_year, arrival_time_vec(i), pop_matrix_male, pop_matrix_female, u(i));
    
    time_pop_sim(i) = pop_sim(0);
    delta_pop_sim(i) = pop_sim(1);
    
  }
  return arma::join_rows(time_pop_sim, delta_pop_sim);
}
