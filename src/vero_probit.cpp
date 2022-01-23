// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "RcppArmadillo.h"
using namespace Rcpp;

// Definicao das funcoes de verossimilhanca e prioris
// O modelo é dado por
//
// y_i|f_i = beta f_i + epsilon_i
// f_i ~ N(0, I_q)
//
// y_i é um vetor de dimensao P
// beta é uma matriz de dimensao P x Q
// f_i é um vetor de dimnesao Q
// epsilon_i é um vetor de dimensao P
// i varia de 1 ate N.
// N é o numero de observacoes
// Q é o numero de fatores
// P é o numero de variáveis
//
// Observamos os dados atraves de y_i^*, que é uma categorizacao de y_i.
// y_ij^* = k se alpha_{k-1} <= y_ij <= alpha_k
// O numero de categorias é prefixado em K. Existem K + 1 pontos de corte,
// onde usualmente os pontos de corte 0 e K são tomados como -infty e +infty
//
// j varia entre 1 e P (quantidade de variaveis)
// i varia entre 1 e N (quantidade de observacoes)
// k varia entre 1 e K (quantidade de categorias)
// l varia entre 1 e Q (quantidade de fatores)

//' P(Y_ij <= k) - Probabilidade da i-ésima observação da j-ésima variável estar abaixo da categoria k. (PROBIT)
//'
//' @param alpha   [1 x 1] escalar com o ponto de corte associado a k-ésima categoria.
//' @param beta   [1 x Q] vetor linha com os betas que associam a j-esima variavel e cada um dos Q fatores.
//' @param f      [Q x 1] vetor coluna com os Q fatores associados a i-ésima observacao.
//' @param sigma2 [1 x 1] escalar com a variância associada à j-ésima variável.
//'
//' @return       [1 x 1] ecalar com a probabilidade de y_ij cair abaixo da categoria k.
//'
//' @export
// [[Rcpp::export]]
double delta_ijk_probit(double alpha, arma::rowvec beta, arma::colvec f, double sigma2){
  return R::pnorm((alpha - arma::dot(beta,f))/std::sqrt(sigma2),0.0,1.0,1,0);
}

// [[Rcpp::export]]
double delta_ijk_marginal_probit(double alpha, arma::rowvec beta, double sigma2){
  int q = beta.size();
  double aux_sum_beta = 0;
  for (int l = 0; l <= q - 1; l++){
    aux_sum_beta += pow(beta(l), 2);
  }
  return R::pnorm(alpha/std::sqrt(aux_sum_beta + sigma2),0.0,1.0,1,0);
}

//' P(Y_ij = k) - Probabilidade da i-ésima observação da j-ésima variável estar abaixo da categoria k. (PROBIT)
//'
//' @param alpha   [1 x 2] vetor com os pontos de corte associados a (k-1)-ésima e (k)-ésima categorias.
//' @param beta   [1 x Q] vetor linha com os betas que associam a j-esima variavel e cada um dos Q fatores.
//' @param f      [Q x 1] vetor coluna com os Q fatores associados a i-ésima observacao.
//' @param sigma2 [1 x 1] escalar com a variância associada à j-ésima variável.
//'
//' @return       [1 x 1] ecalar com a probabilidade de y_ij ser exatamente igual acategoria k.
//'
//' @export
// [[Rcpp::export]]
double p_ijk_probit(arma::rowvec alpha, arma::rowvec beta, arma::colvec f, double sigma2){
  return delta_ijk_probit(alpha(1), beta, f, sigma2) - delta_ijk_probit(alpha(0), beta, f, sigma2);
}

// [[Rcpp::export]]
double p_ijk_marginal_probit(arma::rowvec alpha, arma::rowvec beta, double sigma2){
  return delta_ijk_marginal_probit(alpha(1), beta, sigma2) - delta_ijk_marginal_probit(alpha(0), beta, sigma2);
}

//' P(Y_j = y_j) - (log) Verossimilhanca para a j ésima variavel.
//'
//' @param y      [1 x N] vetor com as n obersavoes para a j-esima variavel.
//' @param alpha   [1 x (K+1)] vetor com os todos os K + 1 pontos de corte.
//' @param beta   [1 x Q] vetor linha com os betas que associam a j-esima variavel e cada um dos Q fatores.
//' @param f      [Q x N] matriz f completa.
//' @param sigma2 [1 x 1] escalar com a variância associada à j-ésima variável.
//'
//' @return       [1 x 1]  ecalar com a log verossimilhanca para todas as observacoes da variavel j.
//'
//' @export
// [[Rcpp::export]]
double l_vero_j_probit(arma::rowvec y,arma::rowvec alpha, arma::rowvec beta, arma::mat f, double sigma2){
  int n = y.size();
  arma::rowvec alpha_aux(2);
  double sum = 0;
  for(int i = 0; i <= n - 1; i++){
    alpha_aux(0) = alpha(y(i) - 1);
    alpha_aux(1) = alpha(y(i));
    sum += std::log(p_ijk_probit(alpha_aux, beta, f.col(i), sigma2));
  }
  return sum;
}

// Essa aqui tá errada. Quando não estou condicionando em f, as variáveis não são independentes
// [[Rcpp::export]]
double l_vero_j_marginal_probit(arma::rowvec y,arma::rowvec alpha, arma::rowvec beta, double sigma2){
  int n = y.size();
  arma::rowvec alpha_aux(2);
  double sum = 0;
  for(int i = 0; i <= n - 1; i++){
    alpha_aux(0) = alpha(y(i) - 1);
    alpha_aux(1) = alpha(y(i));
    sum += std::log(p_ijk_marginal_probit(alpha_aux, beta, sigma2));
  }
  return sum;
}

//' P(Y_i = y_i) - (log) Verossimilhanca para a i-esima observacao.
//'
//' @param y      [P x 1] vetor com as p variaveis para a i-esima observacao.
//' @param alpha   [P x (K+1)]  Matriz alpha com os todos os K + 1 pontos de corte. Pontos de
//' corte podem diferir para diferentes variaveis, embora nao seja usual.
//' @param beta   [P x Q] matriz beta completa.
//' @param f      [Q x 1] vetor coluna com os Q fatores associados a i-ésima observacao.
//' @param sigma2 [P x 1] vetor sigma2 completo.
//'
//' @return       [1 x 1]  ecalar com a log verossimilhanca para todas as variaveis da observacao i.
//'
//' @export
// [[Rcpp::export]]
double l_vero_i_probit(arma::colvec y,arma::mat alpha, arma::mat beta, arma::colvec f, arma::colvec sigma2){
  int p = y.size();
  arma::rowvec alpha_aux(2);
  double sum = 0;
  for(int j = 0; j <= p - 1; j++){
    alpha_aux(0) = alpha(j, y(j) - 1);
    alpha_aux(1) = alpha(j, y(j));
    sum += std::log(p_ijk_probit(alpha_aux, beta.row(j), f, sigma2(j)));
  }
  return sum;
}

//' P(Y = y) - (log) Verossimilhanca (completa).
//'
//' @param y      [P x N] matriz de dados completa.
//' @param alpha   [P x (K+1)] Matriz alpha com os todos os K + 1 pontos de corte. Pontos de
//' corte podem diferir para diferentes variaveis, embora nao seja usual.
//' @param beta   [P x Q] matriz beta completa.
//' @param f      [Q x N] matriz f completa.
//' @param sigma2 [P x 1] vetor sigma2 completo.
//'
//' @return       [1 x 1] ecalar com a log verossimilhanca completa.
//'
//' @export
// [[Rcpp::export]]
double l_vero_probit(arma::mat y,arma::mat alpha, arma::mat beta, arma::mat f, arma::colvec sigma2){
  int p = y.n_rows;
  double sum = 0;
  for(int j = 0; j <= p - 1; j++){
    sum += l_vero_j_probit(y.row(j), alpha.row(j), beta.row(j), f, sigma2(j));
  }
  return sum;
}
// Essa aqui tá errada. Quando não estou condicionando em f, as variáveis não são independentes
// [[Rcpp::export]]
double l_vero_marginal_probit(arma::mat y,arma::mat alpha, arma::mat beta, arma::colvec sigma2){
  int p = y.n_rows;
  double sum = 0;
  for(int j = 0; j <= p - 1; j++){
    sum += l_vero_j_marginal_probit(y.row(j), alpha.row(j), beta.row(j), sigma2(j));
  }
  return sum;
}
