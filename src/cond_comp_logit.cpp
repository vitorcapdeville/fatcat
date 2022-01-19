#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins("cpp11")]]

// Trocar isso aqui por arma::dot
double innerProduct(NumericVector x,
                    NumericVector y) {
  double val = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return val;
}


// Passar essas funcoes para usar o armadillo e fazer o mesmo que foi feito
// no probit

// [[Rcpp::export]]
NumericVector probLogit(NumericVector alfa, NumericVector beta, NumericVector f, double sigma2){
  int k = alfa.size() + 1;
  NumericVector prob(k);
  prob[0] = R::plogis(alfa[0],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
  double sum = prob[0];
  for(int s = 1; s <= k - 2; s++){
    prob[s] = R::plogis(alfa[s],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0)-R::plogis(alfa[s-1],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
    sum += prob[s];
  }
  prob[k-1] = 1 - sum;
  return prob;
}


// [[Rcpp::export]]
double logcondcompbetajLogit(NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2, NumericMatrix y,int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i) - 1;
    aux += std::log(probLogit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = 0;
  int t=0;
  while(t<j-1 && t < beta.size()-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
    t += 1;
  }
  if(t==j-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
  }
  return aux+aux2;
}


// [[Rcpp::export]]
double logcondcompfiLogit(NumericVector f, NumericMatrix beta, NumericMatrix alfa, NumericVector sigma2,NumericMatrix y,int i){
  double aux = 0;
  int j = sigma2.size();
  int p = f.size();
  for(int l = 0; l <= j-1; l++){
    int ind = y(l,i-1) - 1;
    aux += log(probLogit(alfa(l,_),beta(l,_),f,sigma2[l])[ind]);
  }
  double aux2 = 0;
  for(int t=0; t<=p-1;t++){
    aux2 -= std::pow(f[t],2)/2;
  }
  return(aux + aux2);
}


// [[Rcpp::export]]
double logcondcompsigma2jLogit (NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2,NumericMatrix y, int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i)-1;
    aux += log(probLogit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = aux + 0.999*std::log(sigma2) - 0.001/sigma2;
  return(aux2);
}


