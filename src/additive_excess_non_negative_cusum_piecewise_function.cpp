#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "header_rcpp_functions.h"

// Function to calculate the CUSUM chart with the additive alternative, assuming that the excess hazard is non-negative.

arma::vec out_of_control_add_excess_non_neg_cumulative_excess_hazard_piecewise_func(
    arma::vec& partition_t_vec,
    arma::vec& baseline_vec,
    arma::vec& t_vec, 
    arma::vec& linear_predictor_vec,
    double gamma
) {
  
  int n_elements = linear_predictor_vec.size();
  
  arma::vec store_vec(n_elements);
  
  for (int i = 0; i < n_elements; i++) {
    
    store_vec[i] = out_of_control_add_excess_non_neg_cumulative_excess_hazard_piecewise_func_single(partition_t_vec, baseline_vec, t_vec[i], linear_predictor_vec[i], gamma);
    
  }
  
  return store_vec;
  
}

//' Function to calculate the CUSUM chart at different time points given by t_grid
//' 
//' Function to calculate the CUSUM chart at different time points given by t_grid using additive alternative, assuming non-negative excess hazard
//' 
//' @param partition_t_vec A vector that defines the partition of the time period for the piecewise constant baseline hazard.
//' @param baseline_vec A vector of baseline hazard values defined over the time partition. 
//' @param start_year Calendar year of when the monitoring starts. End year is implicitly defined via t_grid. 
//' @param age_vec A vector containing age of the observations when they arrive in the study. 
//' @param gender_vec A vector containing gender indicator - 1 for male, 2 for female. 
//' @param x_matrix A matrix containing the predictor variables (need to contain the same amount of columns as the dimension of beta_vec, the estimated parameter vector).
//' @param time_obs_vec A vector containing the observed time of the observations. 
//' @param arrival_time_vec A vector containing the amount of time units after start year the observations arrive. Should be less than or equal to max(t_grid).
//' @param delta_i_vec A vector containing the event indicator at the end of follow-up. 
//' @param beta_vec A vector containing the estimated parameters from the model related to the predictor variables.
//' @param gamma Additive constant for the out-of-control hazard model. 
//' @param t_grid The time grid to calculate the CUSUM chart. 
//' @param pop_data_male A matrix of population table for male. Should follow the matrix format with first column being the year, second column being age and third column being the hazard values.
//' @param pop_data_female Similar as above, but for female. 
//' @param end_year_table The last calendar year appearing in the population table. 
//' @return final_result_matrix A matrix where the first column is the chart value, the second column counts the number of events.
//' @export
// [[Rcpp::export]]
NumericMatrix cusum_add_excess_non_neg_r_t_piecewise(
    arma::vec partition_t_vec,
    arma::vec baseline_vec,
    double start_year,
    arma::vec age_vec,
    arma::vec gender_vec,
    arma::mat x_matrix,
    arma::vec time_obs_vec,
    arma::vec arrival_time_vec,
    arma::vec delta_i_vec,
    arma::vec beta_vec,
    double gamma,
    arma::vec t_grid,
    arma::mat pop_data_male,
    arma::mat pop_data_female,
    int end_year_table
) {
  
  int length_t_grid = t_grid.size();
  NumericVector R(length_t_grid);
  arma::vec pop_hazard_old; // use to store population hazard in order to not look it up at every time point on t_grid as the population hazard for an individual is only needed at the observed time. Only do this for observations that actually experience event at the end of follow-up.
  arma::vec excess_hazard_old; // same as above
  arma::uvec index_pop_old; // similar as above, store the indices of arrived observations.
  IntegerVector n_events(length_t_grid);
  
  // uvec col_year_index = {0};
  // uvec col_age_index = {1};
  arma::uvec col_haz_index = {2};
  
  for (int i = 0; i < length_t_grid; i++) {
    
    arma::uvec index_vec = find(t_grid(i) >= arrival_time_vec); // find the index of observation that has arrived at this time point of the time grid
    
    int length_of_index = index_vec.size();
    
    if (length_of_index == 0) { // if no one has arrived
      
      double r_t = 0;
      R(i) = r_t;
      
    } else {
      
      // Select all the relevant quantities for the arrived observations
      
      arma::vec new_age_vec = age_vec.elem(index_vec); 
      arma::vec new_gender_vec = gender_vec.elem(index_vec);
      arma::mat new_x_matrix = x_matrix.rows(index_vec);
      arma::vec new_time_obs_vec = time_obs_vec.elem(index_vec);
      arma::vec new_arrival_time_vec = arrival_time_vec.elem(index_vec);
      arma::vec new_delta_i_vec = delta_i_vec.elem(index_vec);
      
      arma::vec at_risk_time_vec = x_i_vec_func(new_time_obs_vec, new_arrival_time_vec, t_grid(i)); // calculate at risk time at the given time point.
      
      arma::uvec index_pop_vec;
      
      // Find indices of observations that have arrived since last time point on the grid that will experience an event at the end.
      // This is for the first term of R(t), where the observation only affects the sum if delta_i(t) = 1. 
      
      if (i == 0) {
        
        index_pop_vec = find(t_grid(i) >= arrival_time_vec && delta_i_vec == 1);
        
      } else {
        
        index_pop_vec = find(t_grid(i - 1) < arrival_time_vec && arrival_time_vec <= t_grid(i) && delta_i_vec == 1);
        
      }
      
      int length_of_index_pop = index_pop_vec.size();
      
      if (length_of_index_pop != 0) {
        // Select all the relevant quantities for these specific observations
        arma::vec pop_hazard_vec(length_of_index_pop);
        arma::vec new_age_vec_pop = age_vec.elem(index_pop_vec);
        arma::vec new_gender_vec_pop = gender_vec.elem(index_pop_vec);
        arma::vec new_time_obs_vec_pop = time_obs_vec.elem(index_pop_vec);
        arma::vec new_arrival_time_vec_pop = arrival_time_vec.elem(index_pop_vec);
        
        arma::vec new_age_floor_vec_pop = floor(new_age_vec_pop + new_time_obs_vec_pop); // age and year are integers in population table
        arma::uvec index_older_than_110 = find(new_age_floor_vec_pop > 110);
        new_age_floor_vec_pop.elem(index_older_than_110).fill(110);
        
        arma::vec new_year_floor_vec_pop = floor(start_year + new_arrival_time_vec_pop + new_time_obs_vec_pop);
        arma::uvec index_year_later_than_end_year_table = find(new_year_floor_vec_pop >= (end_year_table + 1));
        new_year_floor_vec_pop.elem(index_year_later_than_end_year_table).fill(end_year_table);
        
        arma::mat new_x_matrix_pop = x_matrix.rows(index_pop_vec);
        
        for (int j = 0; j < length_of_index_pop; j++) {
          
          if (new_gender_vec_pop(j) == 1) {
            
            arma::uvec index_ind = find(pop_data_male.col(1) == new_age_floor_vec_pop(j) && pop_data_male.col(0) == new_year_floor_vec_pop(j));
            arma::mat pop_hazard_mat = pop_data_male.submat(index_ind, col_haz_index);
            pop_hazard_vec(j) = pop_hazard_mat(0);
            
          } else if (new_gender_vec_pop(j) == 2) {
            
            arma::uvec index_ind = find(pop_data_female.col(1) == new_age_floor_vec_pop(j) && pop_data_female.col(0) == new_year_floor_vec_pop(j));
            arma::mat pop_hazard_mat = pop_data_female.submat(index_ind, col_haz_index);
            pop_hazard_vec(j) = pop_hazard_mat(0);
            
          }
          
        }
        
        pop_hazard_old = join_cols(pop_hazard_old, pop_hazard_vec); // store the obtained hazard values for these observations
        
        arma::vec baseline_haz_vec = baseline_excess_hazard_piecewise_func(partition_t_vec, baseline_vec, new_time_obs_vec_pop); // calculate the estimated baseline hazard for these observations at the observed times
        arma::vec excess_hazard_vec = baseline_haz_vec % (exp(new_x_matrix_pop * beta_vec));
        excess_hazard_old = join_cols(excess_hazard_old, excess_hazard_vec);
        
        index_pop_old = join_cols(index_pop_old, index_pop_vec);
        
        // std::cout << "pop_hazard_old:" << pop_hazard_old.t() << std::endl;
        // std::cout << "baseline_haz_vec:" << baseline_haz_vec.t() << std::endl;
        // std::cout << "excess_hazard_old:" << excess_hazard_old.t() << std::endl;
        // std::cout << "index_pop_old:" << index_pop_old.t() << std::endl;
        
        
      }
      
      if (index_pop_old.size() != 0) {
        
        // First part of R(t)
        arma::vec h1_excess_part_vec = excess_hazard_old + gamma;
        arma::uvec negative_haz_indices = arma::find(h1_excess_part_vec < 0);
        h1_excess_part_vec.elem(negative_haz_indices).fill(0); // assuming that the excess hazard is non-negative
        arma::vec h1_vec = pop_hazard_old + h1_excess_part_vec; // out-of-control hazard for observations experiencing event at the end of follow up
        
        arma::vec h0_vec = pop_hazard_old + excess_hazard_old; // similar to above, but in-control
        
        arma::vec log_haz_ratio_vec = log(h1_vec / h0_vec); // log-likelihood ratio between out-of-control and in-control
        
        arma::vec relevant_time_obs_vec = time_obs_vec.elem(index_pop_old);
        arma::vec relevant_arrival_time_vec = arrival_time_vec.elem(index_pop_old);
        arma::vec relevant_delta_i_vec = delta_i_vec.elem(index_pop_old);
        
        arma::vec event_vec = delta_i_vec_func(relevant_time_obs_vec, relevant_arrival_time_vec, relevant_delta_i_vec, t_grid(i));
        
        // Second part of R(t)
        
        arma::vec linear_predictor_vec = new_x_matrix * beta_vec;
        
        arma::vec H0_vec = cumulative_baseline_excess_hazard_piecewise_func(partition_t_vec, baseline_vec, at_risk_time_vec) % exp(linear_predictor_vec);
        
        arma::vec H1_vec = out_of_control_add_excess_non_neg_cumulative_excess_hazard_piecewise_func(partition_t_vec, baseline_vec, at_risk_time_vec, linear_predictor_vec, gamma);
        
        arma::vec diff_cumulative_haz_vec = H1_vec - H0_vec; 
        
        // Final
        
        n_events(i) = sum(event_vec);
        
        double r_t = sum(event_vec % log_haz_ratio_vec) - sum(diff_cumulative_haz_vec);
        
        R(i) = r_t;
        
      } else {
        
        arma::vec linear_predictor_vec = new_x_matrix * beta_vec;
        
        arma::vec H0_vec = cumulative_baseline_excess_hazard_piecewise_func(partition_t_vec, baseline_vec, at_risk_time_vec) % exp(linear_predictor_vec);
        
        arma::vec H1_vec = out_of_control_add_excess_non_neg_cumulative_excess_hazard_piecewise_func(partition_t_vec, baseline_vec, at_risk_time_vec, linear_predictor_vec, gamma);
        
        arma::vec diff_cumulative_haz_vec = H1_vec - H0_vec; 
        
        n_events(i) = 0;
        
        double r_t = -sum(diff_cumulative_haz_vec);
        
        R(i) = r_t;
        
      }
      
    }
    
  }
  
  NumericVector cumulative_min_R = cummin(R);
  NumericVector chart_vec = R - cumulative_min_R;
  
  NumericMatrix final_result_matrix(length_t_grid, 2);
  
  final_result_matrix(_, 0) = chart_vec;
  final_result_matrix(_, 1) = n_events;
  
  return final_result_matrix;
  
}
