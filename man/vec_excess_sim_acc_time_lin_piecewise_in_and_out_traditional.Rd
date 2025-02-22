% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{vec_excess_sim_acc_time_lin_piecewise_in_and_out_traditional}
\alias{vec_excess_sim_acc_time_lin_piecewise_in_and_out_traditional}
\title{Function to simulate excess times during a time period where the hazard switches from in-control to out-of-control}
\usage{
vec_excess_sim_acc_time_lin_piecewise_in_and_out_traditional(
  partition_t_vec,
  baseline_vec,
  max_follow_up_vec,
  u_vec,
  arrival_time_vec,
  linear_predictor_vec,
  k,
  eta
)
}
\arguments{
\item{partition_t_vec}{A vector that defines the partition of the time period for the piecewise constant baseline hazard.}

\item{baseline_vec}{A vector of baseline hazard values defined over the time partition.}

\item{max_follow_up_vec}{A vector corresponding to the max follow up times of the observations.}

\item{u_vec}{A vector of samples from Unif(0, 1) for inverse transform algorithm.}

\item{arrival_time_vec}{A vector corresponding to the arrival time after the start of monitoring of the observations.}

\item{linear_predictor_vec}{the matrix vector product of the covariate matrix and the regression coefficients.}

\item{k}{Linear time acceleration factor for the out-of-control hazard model.}

\item{eta}{The time point after monitoring start where the hazard switches to the out-of-control.}
}
\value{
time_excess_sim A vector of simulated excess times based on the input of the observations.
}
\description{
Function to simulate excess times during a time period where the hazard switches from in-control to out-of-control using linear accelerated time alternative
}
