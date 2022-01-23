// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins("cpp11")]]

#include "RcppArmadillo.h"
#include "RcppTN.h"
#include <progress.hpp>
#include <eta_progress_bar.hpp>
#include <indica.h>
#include <vero_probit.h>
using namespace Rcpp;

arma::mat passo_beta_probit(arma::mat y, arma::mat alpha, arma::mat beta, arma::mat f, arma::colvec sigma2, double C0, double sdpropbeta, double sdpropbeta2){

  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;

  int p = beta.n_rows;
  int q = beta.n_cols;

  arma::mat ret(p,q);
  arma::rowvec betaprop(q);
  arma::vec auxvec(2);

  for(int j=0; j <= p - 1; j++){
    for(int l=0; l <= q - 1; l++){
      betaprop(l) = R::rnorm(beta(j,l),sdpropbeta)*indica(j>l) + RcppTN::rtn1(beta(j,l),sdpropbeta2,0,R_PosInf)*indica(j==l);
    }
    A1 = l_vero_j_probit(y.row(j), alpha.row(j), betaprop,    f, sigma2(j)) + l_priori_beta_j(betaprop,    j, C0);
    A2 = l_vero_j_probit(y.row(j), alpha.row(j), beta.row(j), f, sigma2(j)) + l_priori_beta_j(beta.row(j), j, C0);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u < pa){
      ret.row(j) = betaprop;
    }else{
      ret.row(j) = beta.row(j);
    }
  }
  return ret;
}

arma::colvec passo_sigma2_probit(arma::mat y, arma::mat alpha, arma::mat beta, arma::mat f, arma::colvec sigma2, double a, double b, double sdpropsigma2){

  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;

  int p = beta.n_rows;

  arma::colvec ret(p);
  double sigma2prop = 0.0;
  arma::vec auxvec(2);

  for(int j=0; j <= p-1; j++){
    sigma2prop = RcppTN::rtn1(sigma2(j),sdpropsigma2,0,R_PosInf);
    A1 = l_vero_j_probit(y.row(j), alpha.row(j), beta.row(j), f, sigma2prop) + l_priori_sigma2_j(sigma2prop, a, b);
    A2 = l_vero_j_probit(y.row(j), alpha.row(j), beta.row(j), f, sigma2(j))  + l_priori_sigma2_j(sigma2(j) , a, b);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u < pa){
      ret(j) = sigma2prop;
    } else {
      ret(j) = sigma2(j);
    }
  }
  return ret;
}

arma::mat passo_f_probit(arma::mat y, arma::mat alpha, arma::mat beta, arma::mat f, arma::colvec sigma2, double sdpropf){
  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;

  int q = beta.n_cols;
  int n = f.n_cols;

  arma::mat ret(q,n);
  arma::colvec fprop(q);
  // arma::colvec aux6(q);
  NumericVector auxvec(2);

  for (int i = 0; i <= n - 1; i++) {
    for (int l = 0; l <= q - 1; l++){
      fprop(l) = R::rnorm(f(l,i), sdpropf);
    }
    // aux6 = Rcpp::rnorm(q)*sdpropf;
    // fprop = f.col(i) + aux6;
    A1 = l_vero_i_probit(y.col(i), alpha, beta, fprop, sigma2)    + l_priori_f_i(fprop);//logcondcompfiProbit(fprop,beta,alpha,sigma2, y,i + 1);
    A2 = l_vero_i_probit(y.col(i), alpha, beta, f.col(i), sigma2) + l_priori_f_i(f.col(i));//logcondcompfiProbit(f.col(i),beta,alpha,sigma2,y,i + 1);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u<pa){
      ret.col(i) = fprop;
    }else{
      ret.col(i) = f.col(i);
    }
  }
  return ret;
}

//' Algoritmo de MH para gerar amostras da posteriori do modelo fatorial categórico ordenado probit
//'
//' @param y      [P x N] matriz de dados completa.
//' @param nit    [1 x 1] numero de interacoes
//' @param beta   [P x Q] valor inicial da matriz de cargas fatoriais.
//' @param f      [Q x N] valor inicial da matriz de fatores.
//' @param sigma2 [P x 1] valor inicial do vetor de variâncias.
//' @param alpha  [P x K + 1] matriz com os valores fixados dos pontos de corte.
//' @param C0     [1 x 1] variância das prioris para beta
//' @param a,b   [1 x 1] parametros de shape e scale da priori de sigma2
//' @param sdpropbeta,sdpropbeta2,sdpropf,sdpropsigma2 desvio padrao das propostas
//' @param display_progress TRUE para exibir o indicador de progresso.
//'
//' @export
// [[Rcpp::export]]
List algoritmo_probit(arma::mat y, int nit , arma::mat beta, arma::mat f, arma::colvec sigma2, arma::mat alpha,
                      double C0 = 10000, double a = 0.001, double b = 0.001,
                      double sdpropbeta = 0.05, double sdpropbeta2 = 0.075, double sdpropf = 0.4, double sdpropsigma2 = 0.05,
                      bool display_progress = true){

  int p = y.n_rows;
  int n = y.n_cols;
  int q = f.n_rows;

  arma::cube chain_beta(p, q, nit);
  arma::cube chain_f(q, n, nit);
  arma::mat chain_sigma2(p, nit);

  chain_beta.slice(1) = beta;
  chain_f.slice(1) = f;
  chain_sigma2.col(1) = sigma2;

  arma::cube chain_export_beta;
  arma::cube chain_export_f;
  arma::mat chain_export_sigma2;

  ETAProgressBar pb;
  Progress pro(nit, display_progress,pb);

  // Progress pro(nit, display_progress);
  for(int it=1; it <= nit-1;it++){

    if (Progress::check_abort()){
      chain_export_beta = chain_beta.slices(1, it - 1);
      chain_export_f = chain_f.slices(1, it - 1);
      chain_export_sigma2 = chain_sigma2.cols(1, it - 1);
      return List::create(_["beta"] = chain_export_beta, _["f"] = chain_export_f, _["sigma2"] = chain_export_sigma2);
    }

    pro.increment();

    chain_beta.slice(it) = passo_beta_probit(y, alpha, chain_beta.slice(it-1), chain_f.slice(it-1), chain_sigma2.col(it-1), C0, sdpropbeta, sdpropbeta2);
    // chain_beta.slice(it) = beta;

    chain_sigma2.col(it) = passo_sigma2_probit(y, alpha, chain_beta.slice(it), chain_f.slice(it-1), chain_sigma2.col(it-1), a, b, sdpropsigma2);
    // chain_sigma2.col(it) = sigma2;

    chain_f.slice(it) = passo_f_probit(y, alpha, chain_beta.slice(it), chain_f.slice(it-1), chain_sigma2.col(it), sdpropf);
    // chain_f.slice(it) = f;
  }
  chain_export_beta = chain_beta.tail_slices(nit - 1);
  chain_export_f = chain_f.tail_slices(nit - 1);
  chain_export_sigma2 = chain_sigma2.tail_cols(nit - 1);
  List res = List::create(_["beta"] = chain_export_beta, _["f"] = chain_export_f, _["sigma2"] = chain_export_sigma2);
  return res;
}
