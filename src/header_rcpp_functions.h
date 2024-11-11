#ifndef MY_RCPP_FILE_H
#define MY_RCPP_FILE_H

// Header file that contains some functions that can be used in several of the script files in order to not redefine it everytime it is needed.

#include <RcppArmadillo.h>
arma::vec delta_i_vec_func(arma::vec&, arma::vec&, arma::vec&, double&);
arma::vec x_i_vec_func(arma::vec&, arma::vec&, double&);
int find_interval_single(double&, arma::vec&);
double cumulative_baseline_excess_hazard_piecewise_func_single(
    arma::vec&, 
    arma::vec&,
    double&
);
double out_of_control_prop_cumulative_baseline_excess_hazard_piecewise_func_single(
    arma::vec&, 
    arma::vec&,
    double,
    double
);
arma::vec baseline_excess_hazard_piecewise_func(
    arma::vec&, 
    arma::vec&,
    arma::vec&
);
double out_of_control_add_excess_non_neg_cumulative_excess_hazard_piecewise_func_single(
    arma::vec&, 
    arma::vec&,
    double,
    double,
    double
);
arma::vec cumulative_baseline_excess_hazard_piecewise_func(
    arma::vec&, 
    arma::vec&,
    arma::vec&
);
double out_of_control_acc_time_lin_cumulative_excess_hazard_add_piecewise_func_single_non_traditional(
    arma::vec&, 
    arma::vec&,
    double,
    double
);
double out_of_control_acc_time_lin_cumulative_excess_hazard_add_piecewise_func_single_traditional(
    arma::vec&, 
    arma::vec&,
    double,
    double
);
  
#endif