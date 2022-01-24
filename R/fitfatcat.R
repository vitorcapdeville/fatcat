#' Função para gerar amostras da posterior, toda em C++.
#'
#' @param y dados. deve ter dimensao jxn
#' @param q numero de fatores a serem usados no ajuste
#' @param nit numero de interacoes
#' @param burnin numero de interacoes iniciais a serem ignoradas
#' @param lag pegar uma a cada `lag` interacoes, após burnin.
#' @param init_beta valor inicial de beta
#' @param init_sigma2 valor inicial de sigma2
#' @param C0 variância das prioris para beta
#' @param a,b parametros de shape e scale da priori de sigma2
#' @param sdpropbeta desvio padrao da proposta
#' @param sdpropbeta2 desvio padrao da proposta
#' @param sdpropf desvio padrao da proposta
#' @param sdpropsigma2 desvio padrao da proposta
#' @param dist distribuição do erro. atualmente, probit ou logit.
#' @param quiet se TRUE, impede o print do indicar de progresso
#' @param alpha opcional, fornecer valores verdadeiros de alpha.
#'
#' @return uma lista com as cadeias a posteriori
#'
#' @import RcppTN
#'
#' @export
fitfatcat <- function(y, q, nit, burnin = 0, lag = 1, init_beta = 0.5, init_sigma2 = 2, C0 = 100000, a = 0.001, b = 0.001, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075, dist = c("probit", "logit"), quiet = F, alpha = NULL) {
  # Definicoes e valores iniciais
  dist <- match.arg(dist)
  if (dist != "probit") stop("Distribuicao logistica em construcao.")

  p <- nrow(y)
  n <- ncol(y)

  beta <- array(0, dim = c(p, q))
  sigma2 <- array(NA, dim = c(p))
  f <- array(NA, dim = c(q, n))
  beta[, 1] <- init_beta
  sigma2 <- rep(init_sigma2, p)
  f <- matrix(0.5, nrow = q, ncol = n)

  if (is.null(alpha)) alpha <- cbind(-Inf,psych::polychoric(t(y))$tau, Inf)

  # Algoritmo
  tempo1 <- Sys.time()
  res <- algoritmo_probit(y, nit + 1, beta, f, sigma2, alpha, C0, a, b, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, !quiet)

  # if (dist == "probit") {
  #   res <- algoritmo_probit(y, nit, beta, sigma2, f, alpha, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, q, p, !quiet)
  # } else {
  #   res <- algoritmo_logit(y, nit, beta, sigma2, f, alpha, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, q, p, !quiet)
  # }

  trim <- trunc(seq(from = burnin + 1, to = min(nit,dim(res$beta)[3]), by = lag))
  res$beta <- res$beta[, , trim, drop = F]
  res$sigma2 <- res$sigma2[, trim, drop = F]
  res$f <- res$f[, , trim, drop = F]

  tempo2 <- Sys.time()

  if (!quiet) message("Executado em ", round(difftime(tempo2, tempo1), 2), " ", units(difftime(tempo2, tempo1)))

  return(res)
}
