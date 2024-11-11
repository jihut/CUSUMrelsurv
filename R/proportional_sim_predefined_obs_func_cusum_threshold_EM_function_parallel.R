`%dopar%` <- foreach::`%dopar%`

# Function to simulate the max values of the CUSUM chart under in-control to estimate the threshold based on this result with the proportional alternative.  
# Here, datasets used to calculate the CUSUM charts are simulated using the true covariate distribution. 

#' Function to simulate threshold for CUSUM with proportional alternative under piecewise constant baseline hazard. 
#' 
#' Function to simulate threshold for CUSUM with proportional alternative under piecewise constant baseline hazard (with covariates simulating via prespecified instructions).
#' @param sim_obs_func A function that simulates covariates etc. under prespecified (true) covariate distribution. 
#' @param base_em_model Fitted model object from relsurv using the semi-parametric model.
#' @param cumulative_baseline_haz_fit A smoothed fit of the estimated cumulative baseline hazard output of the semi-parametric model. 
#' @param baseline_haz_fit A smoothed fit of the estimated baseline hazard output of the semi-parametric model. 
#' @param min_cumulative_baseline_event_time The smallest time to event that was used in the semi-parametric model for estimation.
#' @param n_iterations Number of simulated CUSUM charts under in-control scenario for calculating threshold. 
#' @param n_cores Number of cores to run in parallel.
#' @param rho_vec A vector of values of the proportional constant for the out-of-control hazard model.  
#' @param pop_data_male A matrix of population table for male. Should follow the matrix format with first column being the year, second column being age and third column being the hazard values.
#' @param pop_data_female Similar as above, but for female. 
#' @param end_year_table The last calendar year appearing in the population table. 
#' @param random_state Seed value for reproducibility. 
#' @return store_matrix A matrix with the number of columns equal to the number of rho values. The number of rows is 'n_iterations'. The matrix stores the max value of all the simulated charts, which can be used to estimate the threshold.
#' @export
sim_predefined_obs_func_cusum_prop_threshold_EM_function_parallel <- function(
    sim_obs_func,
    base_em_model = NULL,
    cumulative_baseline_haz_fit,
    baseline_haz_fit,
    min_cumulative_baseline_event_time,
    n_iterations,
    n_cores = NULL,
    rho_vec,
    pop_data_male,
    pop_data_female,
    end_year_table,
    random_state = NULL
) {
  
  # General baseline hazard function based on the smoothed fit based on the estimated baseline and cumulative baseline hazard from the semi-parametric model by Perme, Henderson and Stare (2009)
  
  baseline_haz_general_func <- function(t){ # t in years 
    baseline_value <- predict(baseline_haz_fit, t)$y
    final_baseline_value <- ifelse(baseline_value >= 0, baseline_value, 0)
    final_baseline_value
  }
  
  # Similar for the cumulative baseline
  
  cumulative_baseline_haz_general_func <- function(t){ # t in years 
    baseline_value <- ifelse(t > min_cumulative_baseline_event_time,
                             predict(cumulative_baseline_haz_fit, t)$y,
                             predict(cumulative_baseline_haz_fit, min_cumulative_baseline_event_time)$y / min_cumulative_baseline_event_time * t) # linear interpolation if the evaluated time is less than the smallest time to event in the estimated model. 
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
  
  if (!is.null(random_state)) {
    doRNG::registerDoRNG(seed = random_state)
  }
  
  store_matrix <- foreach::foreach(
    i = 1:n_iterations,
    .combine = "rbind"
  ) %dopar% {
    
    store_vec <- numeric(n_elements_rho_vec)
    
    if (is.null(base_em_model)) {
      sim_obs_run <- sim_obs_func()
    } else {
      sim_obs_run <- sim_obs_func(base_em_model) # relevant in case the estimated regression coefficients are used, then the model object from relsurv needs to be an input too for completeness
    }
    
    n_obs <- sim_obs_run$n_obs
    
    u_vec_excess <- runif(n_obs)
    u_vec_pop <- runif(n_obs)
    
    excess_times <- vec_excess_sim_EM(
      predict_cumulative_baseline_haz_general_func = predict_cumulative_baseline_haz_general_func,
      min_cumulative_baseline_event_time = min_cumulative_baseline_event_time,
      u_vec = u_vec_excess,
      linear_predictor_vec = sim_obs_run$x_matrix %*% sim_obs_run$beta_vec,
      max_follow_up_vec = sim_obs_run$max_follow_up_vec
    )
    
    matrix_time_pop_sim <- vec_pop_sim(
      age_vec = sim_obs_run$age_vec,
      gender_vec = sim_obs_run$gender_vec,
      start_year = sim_obs_run$start_year,
      end_year = sim_obs_run$end_year,
      table_end_year = end_year_table,
      arrival_time_vec = sim_obs_run$arrival_time_vec,
      pop_data_male,
      pop_data_female,
      u_vec_pop
    )
    pop_times <- matrix_time_pop_sim[, 1]
    delta_p_vec <- matrix_time_pop_sim[, 2]
    censoring_times <-
      pmin(rexp(n_obs, rate = sim_obs_run$interim_censoring_rate),
           sim_obs_run$max_follow_up_vec)
    observed_times <- pmin(excess_times, pop_times, censoring_times)
    delta_i_vec <-
      pmax(delta_p_vec, as.numeric(excess_times < pop_times)) * as.numeric(censoring_times >
                                                                             observed_times)
    
    for (j in 1:n_elements_rho_vec) {
      cusum_chart <- cusum_prop_r_t_EM_smoothing_spline(
        cumulative_baseline_haz_vec_general_func = cumulative_baseline_haz_general_func,
        baseline_haz_vec_general_func = baseline_haz_general_func,
        start_year = sim_obs_run$start_year,
        age_vec = sim_obs_run$age_vec,
        gender_vec = sim_obs_run$gender_vec,
        x_matrix = sim_obs_run$x_matrix,
        time_obs_vec = observed_times,
        arrival_time_vec = sim_obs_run$arrival_time_vec,
        delta_i_vec = delta_i_vec,
        beta_vec = sim_obs_run$beta_vec,
        rho = rho_vec[j],
        t_grid = sim_obs_run$t_grid,
        pop_data_male = pop_data_male,
        pop_data_female = pop_data_female,
        end_year_table = end_year_table
      )
      
      store_vec[j] <- max(cusum_chart[, 1])
    }
    
    store_vec
    
  }
  
  parallel::stopCluster(init_cluster)
  
  colnames(store_matrix) <- paste("rho_", rho_vec, sep = "")
  
  store_matrix
  
}