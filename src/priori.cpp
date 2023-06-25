// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "RcppArmadillo.h"
using namespace Rcpp;

// Definicao das prioris.
//
// j varia entre 1 e P (quantidade de variaveis)
// i varia entre 1 e N (quantidade de observacoes)
// k varia entre 1 e K (quantidade de categorias)
// l varia entre 1 e Q (quantidade de fatores)


// (log) Priori para beta_jl
//
// @param beta [1 x 1] valor do beta_jl
// @param j    [1 x 1] escalar indicando qual é o j associado a beta_jl
// @param l    [1 x 1] escalar indicando qual é o l associado a beta_jl
// @param C0   [1 x 1] variancia da priori de beta.
//
// @details
// As prioris para beta são:
// beta_jl ~ N(0, C0), se j < l
// beta_jl ~ N(0, C0)I(0, +infty), se j == l
// beta_jl = 0 se j > l
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo) de beta_jl
//

double l_priori_beta_jl(double beta, int j, int l, double C0 = 10000){
  double aux;
  if (l <= j) {
    // Nucleo da normal e da normal truncada no 0 é o mesmo (truncar no 0 = dividir por 0.5, se a media e 0)
    // O certo aqui era ter um outro if pra l == j e dizer que se l == j e beta < 0, entao a densidade
    // é zero, mas, pra ser mais eficiente, eu posso não fazer esse if e garantir que se l == j eu não vou
    // propor beta negativo.
    aux =  (-1/(2 * C0))*std::pow(beta, 2);
  } else {
    aux = 0.0; // Se j > l, beta_jl = 0 com prob 1. log(1) = 0
  }
  return aux;
}

// (log) Priori para beta_j.
//
// @param beta [1 x Q] vetor linha com os betas que associam a j-esima variavel e cada um dos Q fatores.
// @param j    [1 x 1] escalar indicando qual é o j associado ao beta_j fornecido.
// @param C0   [1 x 1] variancia da priori de beta.
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_beta_j(arma::rowvec beta, int j, double C0 = 10000){
  double sum = 0;
  int q = beta.n_elem;
  for(int l = 1; l <= q; l++) {
    sum += l_priori_beta_jl(beta(l - 1), j, l, C0);
  }
  return sum;
}

// (log) Priori completa de beta.
//
// @param beta [P x Q] matriz beta completa.
// @param C0   [1 x 1] variancia da priori de beta.
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_beta(arma::mat beta, double C0 = 10000){
  double sum = 0;
  int p = beta.n_rows;
  for(int j = 1; j <= p; j++) {
    sum += l_priori_beta_j(beta.row(j - 1), j, C0);
  }
  return sum;
}

// (log) Priori para f_li
//
// @param f [1 x 1] valor de f_li
//
// @details
// A priori para f é
// f_li ~ N(0, 1)
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_f_li(double f){
  double aux = -0.5*std::pow(f, 2);
  return aux;
  // return R::dnorm(f, 0.0, 1.0, true);
}

// (log) Priori para f_i
//
// @param f [Q x 1] vetor coluna com os Q fatores associados a i-ésima observacao.
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_f_i(arma::colvec f){
  double sum = 0;
  int q = f.n_elem;
  for(int l = 1; l <= q; l++) {
    sum += l_priori_f_li(f(l - 1));
  }
  return sum;
}

// (log) Priori para f
//
// @param f [Q x N] matriz f completa.
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_f(arma::mat f){
  int n = f.n_cols;
  double sum = 0;
  for(int i = 1; i <= n; i++) {
    sum += l_priori_f_i(f.col(i - 1));
  }
  return sum;
}


// (log) Priori para sigma2_j
//
// @param sigma2 [1 x 1] valor da variancia da j-ésima variavel
// @param a      [1 x 1] parametro de shape da IG.
// @param b      [1 x 1] parametro de scale da IG.
//
// @details
// A priori para sigma2_j é
// sigma2_j ~ IG(a, b)
// f(sigma2_j) propto sigma2_j^(- a - 1) exp (-b/sigma2_j)
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_sigma2_j (double sigma2, double a = 0.001, double b = 0.001){
  // Essa priori deveria ter um if dizendo q se sigma2 < 0 a densidade e zero, mas isso
  // pode ser evitado se nunca sera proposto sigma2 < 0.
  double aux = (- a - 1) * std::log(sigma2) - b/sigma2;
  return aux;
}

// (log) Priori para sigma2
//
// @param sigma2 [P x 1] vetor com as variancias das p variáveis.
// @param a      [1 x 1] parametro de shape da IG.
// @param b      [1 x 1] parametro de scale da IG.
//
// @return [1 x 1] escalar com o valor do log da densidade (nucleo)
//

double l_priori_sigma2 (arma::colvec sigma2, double a = 0.001, double b = 0.001){
  int p = sigma2.n_elem;
  double sum = 0;
  for(int j = 1; j <= p; j++) {
    sum += l_priori_sigma2_j(sigma2(j - 1), a, b);
  }
  return sum;
}
