#' Função para gerar amostras da posterior, toda em C++.
#'
#' @param y dados. deve ter dimensao jxn
#' @param p numero de fatores a serem usados no ajuste
#' @param nit numero de interacoes
#' @param burin numero de interacoes iniciais a serem ignoradas
#' @param lag pegar uma a cada `lag` interacoes, após burnin.
#' @param init_beta valor inicial de beta
#' @param init_sigma2 valor inicial de sigma2
#' @param sdpropbeta desvio padrao da proposta
#' @param sdpropbeta2 desvio padrao da proposta
#' @param sdpropf desvio padrao da proposta
#' @param sdpropsigma2 desvio padrao da proposta
#' @param dist distribuição do erro. atualmente, probit ou logit.
#'
#' @return uma lista com as cadeias a posteriori
#'
#' @import RcppTN
#'
#' @export
fitfatcat <- function(y, p, nit, burnin = nit * 0.05, lag = 1, init_beta = 0.5, init_sigma2 = 2, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075, dist = c("probit", "logit")) {
  # Definicoes e valores iniciais
  dist <- match.arg(dist)
  j <- nrow(y)
  n <- ncol(y)
  beta <- array(0, dim = c(j, p, nit))
  sigma2 <- array(NA, dim = c(j, nit))
  f <- array(NA, dim = c(p, n, nit))
  beta[1, 1, 1] <- init_beta
  sigma2[, 1] <- rep(init_sigma2, j)
  f[, , 1] <- matrix(0, nrow = p, ncol = n)
  alfa <- psych::polychoric(t(y))$tau

  # Algoritmo
  tempo1 <- Sys.time()
  if (dist == "probit") {
    res <- algoritmo_probit(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j)
  } else {
    res <- algoritmo_logit(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j)
  }

  trim <- trunc(seq(from = burnin + 1, to = nit, by = lag))

  res$beta <- array(sapply(trim, function(y) lastInd(res$beta, y)), dim = c(j, p, length(trim)))
  res$sigma2 <- array(sapply(trim, function(y) lastInd(res$sigma2, y)), dim = c(j, length(trim)))
  res$f <- array(sapply(trim, function(y) lastInd(res$f, y)), dim = c(p, n, length(trim)))

  tempo2 <- Sys.time()

  message("Executado em ", round(difftime(tempo2, tempo1), 2), " ", units(difftime(tempo2, tempo1)))

  return(res)
}

#' https://stackoverflow.com/questions/44099454/in-r-how-to-index-only-the-last-dimension-of-an-array-when-you-dont-know-how-m
lastInd <- function(x, n) {
  nd <- length(dim(x))
  # uncomment the following line if x could be a plain vector without dim
  # if(nd == 0) nd = 1
  inds <- rep(alist(, )[1], nd)
  inds[nd] <- n
  do.call(`[`, c(list(x), inds))
}
