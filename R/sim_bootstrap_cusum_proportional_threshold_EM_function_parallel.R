`%dopar%` <- foreach::`%dopar%`

# Function to simulate the max values of the CUSUM chart under in-control to estimate the threshold based on this result with the proportional alternative.  
# Here, datasets used to calculate the CUSUM charts are simulated by bootstrapping the dataset used to estimate the in-control model. 

#' Function to simulate threshold for CUSUM with proportional alternative under the semi-parametric (EM) excess hazard model
#' 
#' Function to simulate threshold for CUSUM with proportional alternative under the semi-parametric (EM) excess hazard model (with covariates simulating via bootstrapping)
#' @param cumulative_baseline_haz_fit AA smoothed fit of the estimated cumulative baseline hazard output of the semi-parametric model. 
#' @param baseline_haz_fit A smoothed fit of the estimated baseline hazard output of the semi-parametric model. 
#' @param min_cumulative_baseline_event_time The smallest time to event that was used in the semi-parametric model for estimation.
#' @param n_iterations Number of simulated CUSUM charts under in-control scenario for calculating threshold. 
#' @param n_cores Number of cores to run in parallel.
#' @param start_year Calendar year of when the monitoring starts. End year is implicitly defined via t_grid. 
#' @param estimated_rate_censoring Estimated censoring rate when simulating interim censoring times for observations. 
#' @param age_vec A vector containing age of the observations when they arrive in the study. 
#' @param gender_vec A vector containing gender indicator - 0 for male, 1 for female. 
#' @param x_matrix A matrix containing the predictor variables (need to contain the same amount of columns as the dimension of beta_vec, the estimated parameter vector).
#' @param time_obs_vec A vector containing the observed time of the observations. 
#' @param arrival_time_vec A vector containing the amount of time units after start year the observations arrive. Should be less than or equal to max(t_grid).
#' @param beta_vec A vector containing the estimated parameters from the model related to the predictor variables.
#' @param rho_vec A vector of values of the proportional constant for the out-of-control hazard model.  
#' @param t_grid The time grid to calculate the CUSUM chart. 
#' @param pop_data_male A matrix of population table for male. Should follow the matrix format with first column being the year, second column being age and third column being the hazard values.
#' @param pop_data_female Similar as above, but for female. 
#' @param end_year_table The last calendar year appearing in the population table. 
#' @param random_state Seed value for reproducibility. 
#' @return store_matrix A matrix with the number of columns equal to the number of rho values. The number of rows is 'n_iterations'. The matrix stores the max value of all the simulated charts, which can be used to estimate the threshold.
#' @export
sim_bootstrap_cusum_prop_threshold_EM_function_parallel <- function(
    cumulative_baseline_haz_fit,
    baseline_haz_fit,
    min_cumulative_baseline_event_time,
    n_iterations,
    n_cores = NULL,
    start_year,
    estimated_rate_censoring,
    age_vec,
    gender_vec,
    x_matrix,
    arrival_time_vec,
    beta_vec,
    rho_vec,
    t_grid,
    pop_data_male,
    pop_data_female,
    end_year_table,
    random_state = NULL
) {
  
  baseline_haz_general_func <- function(t){ # t in years 
    baseline_value <- predict(baseline_haz_fit, t)$y
    final_baseline_value <- ifelse(baseline_value >= 0, baseline_value, 0)
    final_baseline_value
  }
  
  cumulative_baseline_haz_general_func <- function(t){ # t in years 
    baseline_value <- ifelse(t > min_cumulative_baseline_event_time,
                             predict(cumulative_baseline_haz_fit, t)$y,
                             predict(cumulative_baseline_haz_fit, min_cumulative_baseline_event_time)$y / min_cumulative_baseline_event_time * t)
    baseline_value
  }
  
  predict_cumulative_baseline_haz_general_func <- function(t){ # t in years 
    predict(cumulative_baseline_haz_fit, t)$y
  }
  
  if (is.null(n_cores)) {
    
    n_cores <- parallel::detectCores() - 2
    
  }
  
  init_cluster <- parallel::makeCluster(n_cores)
  
  doParallel::registerDoParallel(cl = init_cluster)
  
  print(foreach::getDoParRegistered())
  
  n_elements_rho_vec <- length(rho_vec)
  n_observations <- nrow(x_matrix)
  num_of_years <- t_grid[length(t_grid)]
  end_year <- start_year + num_of_years
  
  if (!is.null(random_state)) {
    doRNG::registerDoRNG(seed = random_state)
  }
  
  store_matrix <- foreach::foreach(
    i = 1:n_iterations,
    .combine = "rbind"
  ) %dopar% {
    
    store_vec <- numeric(n_elements_rho_vec)
    
    index_bootstrap <- sample(1:n_observations, replace = T)
    new_arrival_time_vec <- arrival_time_vec[index_bootstrap]
    final_index <- index_bootstrap[which(new_arrival_time_vec < num_of_years)]
    
    new_x_matrix <- x_matrix[final_index, ]
    new_age_vec <- age_vec[final_index]
    new_gender_vec <- gender_vec[final_index]
    new_arrival_time_vec <- arrival_time_vec[final_index]
    new_max_follow_up_vec <-
      num_of_years - new_arrival_time_vec
    
    u_vec_excess <- runif(length(final_index))
    u_vec_pop <- runif(length(final_index))
    
    excess_times <- vec_excess_sim_EM(
      predict_cumulative_baseline_haz_general_func = predict_cumulative_baseline_haz_general_func,
      min_cumulative_baseline_event_time = min_cumulative_baseline_event_time,
      u_vec = u_vec_excess,
      linear_predictor_vec = new_x_matrix %*% beta_vec,
      max_follow_up_vec = new_max_follow_up_vec
    )
    
    matrix_time_pop_sim <- vec_pop_sim(
      new_age_vec,
      new_gender_vec,
      start_year,
      end_year,
      end_year_table,
      new_arrival_time_vec,
      pop_data_male,
      pop_data_female,
      u_vec_pop
    )
    pop_times <- matrix_time_pop_sim[, 1]
    delta_p_vec <- matrix_time_pop_sim[, 2]
    if (estimated_rate_censoring == 0) {
      censoring_times <- new_max_follow_up_vec
    } else {
      censoring_times <-
        pmin(rexp(length(final_index), rate = estimated_rate_censoring),
             new_max_follow_up_vec)
    }
    observed_times <- pmin(excess_times, pop_times, censoring_times)
    delta_i_vec <-
      pmax(delta_p_vec, as.numeric(excess_times < pop_times)) * as.numeric(censoring_times >
                                                                             observed_times)
    
    for (j in 1:n_elements_rho_vec) {
      cusum.chart <- cusum_prop_r_t_EM_smoothing_spline(
        cumulative_baseline_haz_vec_general_func = cumulative_baseline_haz_general_func,
        baseline_haz_vec_general_func = baseline_haz_general_func,
        start_year = start_year,
        age_vec = new_age_vec,
        gender_vec = new_gender_vec,
        x_matrix = new_x_matrix,
        time_obs_vec = observed_times,
        arrival_time_vec = new_arrival_time_vec,
        delta_i_vec = delta_i_vec,
        beta_vec = beta_vec,
        rho = rho_vec[j],
        t_grid = t_grid,
        pop_data_male = pop_data_male,
        pop_data_female = pop_data_female,
        end_year_table = end_year_table
      )
      
      store_vec[j] <- max(cusum.chart[, 1])
    }
    
    store_vec
    
  }
  
  parallel::stopCluster(init_cluster)
  
  colnames(store_matrix) <- paste("rho_", rho_vec, sep = "")
  
  store_matrix
  
}
