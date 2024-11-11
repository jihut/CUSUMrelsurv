`%dopar%` <- foreach::`%dopar%`

# Function to simulate the max values of the CUSUM chart under in-control to estimate the threshold based on this result with the linear accelerated time alternative. 
# Here, datasets used to calculate the CUSUM charts are simulated using the true covariate distribution. 


#' Function to simulate threshold for CUSUM with linear accelerated time alternative under piecewise constant baseline hazard. 
#' 
#' Function to simulate threshold for CUSUM with linear accelerated time alternative under piecewise constant baseline hazard (with covariates simulating via prespecified instructions).
#' @param sim_obs_func A function that simulates covariates etc. under prespecified (true) covariate distribution. 
#' @param n_iterations Number of simulated CUSUM charts under in-control scenario for calculating threshold. 
#' @param n_cores Number of cores to run in parallel.
#' @param k_vec A vector of values of the linear time acceleration factor for the out-of-control hazard model.  
#' @param pop_data_male A matrix of population table for male. Should follow the matrix format with first column being the year, second column being age and third column being the hazard values.
#' @param pop_data_female Similar as above, but for female. 
#' @param end_year_table The last calendar year appearing in the population table. 
#' @param random_state Seed value for reproducibility. 
#' @return store_matrix A matrix with the number of columns equal to the number of rho values. The number of rows is 'n_iterations'. The matrix stores the max value of all the simulated charts, which can be used to estimate the threshold.
#' @export
sim_predefined_obs_func_cusum_acc_time_lin_threshold_piecewise_function_parallel_traditional <- function(
    sim_obs_func,
    n_iterations,
    n_cores = NULL,
    k_vec,
    pop_data_male,
    pop_data_female,
    end_year_table,
    random_state = NULL
) {
  
  if (is.null(n_cores)) {
    
    n_cores <- parallel::detectCores() - 2
    
  }
  
  init_cluster <- parallel::makeCluster(n_cores)
  
  doParallel::registerDoParallel(cl = init_cluster)
  
  print(foreach::getDoParRegistered())
  
  n_elements_k_vec <- length(k_vec)
  
  if (!is.null(random_state)) {
    doRNG::registerDoRNG(seed = random_state)
  }
  
  store_matrix <- foreach::foreach(
    i = 1:n_iterations,
    .combine = "rbind"
  ) %dopar% {
    
    store_vec <- numeric(n_elements_k_vec)
    
    sim_obs_run <- sim_obs_func() # simulate covariates etc. based on the provided function
    
    n_obs <- sim_obs_run$n_obs
    
    u_vec_excess <- runif(n_obs)
    u_vec_pop <- runif(n_obs)
    
    excess_times <- vec_excess_sim_piecewise( # simulate excess times under the (estimated) in-control hazard
      partition_t_vec = sim_obs_run$partition_t_vec,
      baseline_vec = sim_obs_run$baseline_vec,
      max_follow_up_vec = sim_obs_run$max_follow_up_vec,
      u_vec = u_vec_excess,
      linear_predictor_vec = sim_obs_run$x_matrix %*% sim_obs_run$beta_vec
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
    
    for (j in 1:n_elements_k_vec) {
      cusum_chart <- cusum_acc_time_lin_r_t_piecewise_traditional(
        partition_t_vec = sim_obs_run$partition_t_vec,
        baseline_vec = sim_obs_run$baseline_vec,
        start_year = sim_obs_run$start_year,
        age_vec = sim_obs_run$age_vec,
        gender_vec = sim_obs_run$gender_vec,
        x_matrix = sim_obs_run$x_matrix,
        time_obs_vec = observed_times,
        arrival_time_vec = sim_obs_run$arrival_time_vec,
        delta_i_vec = delta_i_vec,
        beta_vec = sim_obs_run$beta_vec,
        k = k_vec[j],
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
  
  colnames(store_matrix) <- paste("k_", k_vec, sep = "")
  
  store_matrix
  
}