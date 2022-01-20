// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppProgress)]]

#include "RcppArmadillo.h"
#include "RcppTN.h"
#include <progress.hpp>
#include <eta_progress_bar.hpp>
#include <indica.h>
using namespace Rcpp;

arma::rowvec probLogit(arma::rowvec alfa, arma::rowvec beta, arma::colvec f, double sigma2){
  int k = alfa.size() + 1;
  NumericVector prob(k);
  prob(0) = R::plogis(alfa(0), arma::dot(beta,f), std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
  double sum = prob(0);
  for(int s = 1; s <= k - 2; s++){
    prob(s) = R::plogis(alfa(s), arma::dot(beta,f), std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0)-R::plogis(alfa(s-1), arma::dot(beta,f), std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
    sum += prob(s);
  }
  prob(k-1) = 1 - sum;
  return prob;
}

double logcondcompbetajLogit(arma::mat f,arma::rowvec beta,arma::rowvec alfa,double sigma2, arma::mat y,int j){
  double aux = 0;
  int n = f.n_cols;
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i) - 1;
    aux += std::log(probLogit(alfa,beta,f.col(i),sigma2)(ind));
  }
  double aux2 = 0;
  int t=0;
  int beta_size = beta.n_elem;
  while(t < j - 1 && t < beta_size - 1){
    aux2 -= (1/2000000)*std::pow(beta(t),2);
    t += 1;
  }
  if(t==j-1){
    aux2 -= (1/2000000)*std::pow(beta(t),2);
  }
  return aux+aux2;
}

double logcondcompfiLogit(arma::colvec f, arma::mat beta, arma::mat alfa, arma::colvec sigma2, arma::mat y,int i){
  double aux = 0;
  int j = sigma2.size();
  int p = f.size();
  for(int l = 0; l <= j-1; l++){
    int ind = y(l,i-1) - 1;
    aux += log(probLogit(alfa.row(l),beta.row(l),f,sigma2(l))(ind));
  }
  double aux2 = 0;
  for(int t=0; t<=p-1;t++){
    aux2 -= std::pow(f(t),2)/2;
  }
  return(aux + aux2);
}

double logcondcompsigma2jLogit (arma::mat f,arma::rowvec beta,arma::rowvec alfa, double sigma2, arma::mat y, int j){
  double aux = 0;
  int n = f.n_cols;
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i)-1;
    aux += log(probLogit(alfa,beta,f.col(i),sigma2)(ind));
  }
  double aux2 = aux + 0.999*std::log(sigma2) - 0.001/sigma2;
  return(aux2);
}

arma::mat passo_beta_logit(int it, int j, int p, arma::mat aux4, double sdpropbeta, double sdpropbeta2, arma::mat f, arma::mat alfa, arma::colvec sigma2, arma::mat y){
  arma::mat auxlb(j,p);
  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;
  arma::rowvec betaprop(p);

  NumericVector auxvec(2);

  for(int l=0; l <= j - 1; l++){
    for(int t=0; t <= p - 1; t++){
      betaprop(t) = R::rnorm(aux4(l,t),sdpropbeta)*indica(l>t) + RcppTN::rtn1(aux4(l,t),sdpropbeta2,0,R_PosInf)*indica(l==t);
    }
    A1 = logcondcompbetajLogit(f,betaprop,alfa.row(l),sigma2(l),y,l+1);
    A2 = logcondcompbetajLogit(f,aux4.row(l),alfa.row(l),sigma2(l),y,l+1);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u<pa){
      auxlb.row(l) = betaprop;
      aux4.row(l) = betaprop;
    }else{
      auxlb.row(l) = aux4.row(l);
    }
  }
  return auxlb;
}

arma::colvec passo_sigma2_logit(int it, int j, arma::colvec aux3, double sdpropsigma2, arma::mat beta, arma::mat f, arma::mat alfa,arma::mat y){
  double sigma2prop = 0.0;
  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;
  NumericVector auxvec(2);

  for(int l=0; l <= j-1; l++){
    sigma2prop = RcppTN::rtn1(aux3(l),sdpropsigma2,0,R_PosInf);
    A1 = logcondcompsigma2jLogit(f,beta.row(l),alfa.row(l),sigma2prop,y,l+1);
    A2 = logcondcompsigma2jLogit(f,beta.row(l),alfa.row(l),aux3(l),y,l+1);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u<pa){
      aux3(l) = sigma2prop;
    }
  }
  return aux3;
}

arma::mat passo_f_logit(int it, int n, int p, arma::mat aux5, double sdpropf, arma::mat beta, arma::colvec sigma2, arma::mat alfa, arma::mat y){
  arma::colvec aux6(p);
  arma::colvec fprop;
  double A1 = 0.0;
  double A2 = 0.0;
  double A = 0.0;
  double pa = 0.0;
  double u = 0.0;
  NumericVector auxvec(2);
  arma::mat auxlf(p,n);

  for(int i = 0; i <= n-1; i++){
    aux6 = Rcpp::rnorm(p)*sdpropf;
    fprop = aux5.col(i) + aux6;
    A1 = logcondcompfiLogit(fprop,beta,alfa,sigma2, y,i + 1);
    A2 = logcondcompfiLogit(aux5.col(i),beta,alfa,sigma2,y,i + 1);
    A  = std::exp(A1 - A2);
    auxvec = {1.0,A};
    pa = min(auxvec);
    u = R::runif(0.0,1.0);
    if(u<pa){
      auxlf.col(i) = fprop;
      aux5.col(i) = fprop;
    }else{
      auxlf.col(i) = aux5.col(i);
    }
  }
  return auxlf;
}

// [[Rcpp::export]]
List algoritmo_logit(arma::mat y, int nit, arma::cube beta, arma::mat sigma2, arma::cube f, arma::mat alfa, double sdpropbeta, double sdpropbeta2, double sdpropf, double sdpropsigma2, int n, int p, int j,bool display_progress = true){

  ETAProgressBar pb;
  Progress pro(nit, display_progress,pb);

  // Progress pro(nit, display_progress);
  for(int it=1; it <= nit-1;it++){

    if (Progress::check_abort())
      return List::create(_["beta"] = beta , _["f"] = f, _["sigma2"] = sigma2);

    pro.increment();

    beta.slice(it) = passo_beta_logit(it, j, p, beta.slice(it-1), sdpropbeta, sdpropbeta2, f.slice(it-1), alfa, sigma2.col(it-1), y);

    sigma2.col(it) = passo_sigma2_logit(it, j, sigma2.col(it-1), sdpropsigma2 , beta.slice(it), f.slice(it-1), alfa, y);

    f.slice(it) = passo_f_logit(it, n, p, f.slice(it-1), sdpropf, beta.slice(it), sigma2.col(it), alfa, y);
  }
  List res = List::create(_["beta"] = beta , _["f"] = f, _["sigma2"] = sigma2);
  return res;
}
