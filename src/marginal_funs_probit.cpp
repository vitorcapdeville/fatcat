// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins("cpp11")]]
//
// #include "RcppArmadillo.h"
// using namespace Rcpp;
//
// arma::rowvec marginalProbProbit(arma::rowvec alfa, arma::rowvec beta, double sigma2){
//   int k = alfa.size() + 1;
//   int j = beta.size();
//   arma::rowvec prob(k);
//   double aux_sum_beta = 0;
//   for (int l = 0; l <= j - 1; l++){
//     aux_sum_beta += pow(beta(l), 2);
//   }
//   prob(0) = R::pnorm(alfa(0)/(std::sqrt(aux_sum_beta + sigma2)),0.0,1.0,1,0);
//   double sum = prob(0);
//   for(int s = 1; s <= k - 2; s++){
//     prob(s) = R::pnorm(alfa(s)/(std::sqrt(aux_sum_beta + sigma2)),0.0,1.0,1,0) -
//       R::pnorm(alfa(s-1)/(std::sqrt(aux_sum_beta + sigma2)),0.0,1.0,1,0);
//     sum += prob(s);
//   }
//   prob(k-1) = 1 - sum;
//   return prob;
// }
//
// //' @export
// // [[Rcpp::export]]
// double marginalloglikjProbit (arma::rowvec beta,arma::rowvec alfa,double sigma2, arma::mat y,int j){
//   double aux = 0;
//   int n = y.n_cols;
//   for(int i = 0; i <= n-1; i++){
//     int ind = y(j-1,i) - 1;
//     aux += std::log(marginalProbProbit(alfa, beta, sigma2)(ind));
//   }
//   return aux;
// }
