#ifndef VERO_PROBIT_H
#define VERO_PROBIT_H

#include "RcppArmadillo.h"

double l_vero_i_probit(arma::colvec y,arma::mat alpha, arma::mat beta, arma::colvec f, arma::colvec sigma2);
double l_vero_j_probit(arma::rowvec y,arma::rowvec alpha, arma::rowvec beta, arma::mat f, double sigma2);
double l_vero_probit(arma::mat y,arma::mat alpha, arma::mat beta, arma::mat f, arma::colvec sigma2);
double l_priori_beta_j(arma::rowvec beta, int j, double C0 = 10000);
double l_priori_beta(arma::mat beta, double C0 = 10000);
double l_priori_sigma2_j (double sigma2, double a = 0.001, double b = 0.001);
double l_priori_f_i(arma::colvec f);
double l_priori_f(arma::mat f);

#endif
