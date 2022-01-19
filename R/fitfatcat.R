#' Função para gerar amostras da posterior, toda em C++.
#'
#' @param y dados. deve ter dimensao jxn
#' @param p numero de fatores a serem usados no ajuste
#' @param nit numero de interacoes
#' @param init_beta valor inicial de beta
#' @param init_sigma2 valor inicial de sigma2
#' @param sdpropbeta desvio padrao da proposta
#' @param sdpropbeta2 desvio padrao da proposta
#' @param sdpropf desvio padrao da proposta
#' @param sdpropsigma2 desvio padrao da proposta
#'
#' @return uma lista com as cadeias a posteriori
#'
#' @import RcppTN
#'
#' @export
fitfatcat <- function(y, p, nit, init_beta = 0.5, init_sigma2 = 2, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075) {
  # Definicoes e valores iniciais
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
  # Incluir barra de progresso e checkforinterrupt
  tempo1 <- Sys.time()

  res <- algoritmo_probit(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j)

  tempo2 <- Sys.time()

  cat("\nExecutado em", round(difftime(tempo2, tempo1),2), units(difftime(tempo2, tempo1)), "\n")

  return(res)
}



#' Função para gerar amostras da posterior, usando condicionais completas em C++ e
#' o algoritmo ainda no R. Mantido aqui para criterios de comparaçao de resultado
#' e performance
#'
#' @param y dados. deve ter dimensao jxn
#' @param p numero de fatores a serem usados no ajuste
#' @param nit numero de interacoes
#' @param init_beta valor inicial de beta
#' @param init_sigma2 valor inicial de sigma2
#' @param sdpropbeta desvio padrao da proposta
#' @param sdpropbeta2 desvio padrao da proposta
#' @param sdpropf desvio padrao da proposta
#' @param sdpropsigma2 desvio padrao da proposta
#'
#' @return uma lista com as cadeias a posteriori
#'
#'
#' @export
fitfatcat_R <- function(y, p, nit, init_beta = 0.5, init_sigma2 = 2, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075) {
  j <- nrow(y)
  n <- ncol(y)
  beta <- array(0, dim = c(j, p, nit))
  sigma2 <- array(NA, dim = c(j, nit))
  f <- array(NA, dim = c(p, n, nit))
  beta[1, 1, 1] <- init_beta
  sigma2[, 1] <- rep(init_sigma2, j)
  f[, , 1] <- matrix(0, nrow = p, ncol = n)
  alfa <- psych::polychoric(t(y))$tau
  tempo1 <- Sys.time()
  if (p == 1) {
    res <- fatcat_p_igual_a_1(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j)
  } else {
    res <- fatcat_p_maior_q_1(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j)
  }
  tempo2 <- Sys.time()
  cat("Executado em ", round(difftime(tempo2, tempo1),2), units(difftime(tempo2, tempo1)))
  return(res)
}

fatcat_p_maior_q_1 <- function(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j) {
  pb <- utils::txtProgressBar(min = 2, max = nit, style = 3)
  for (it in 2:nit) {
    aux4 <- beta[, , it - 1]
    aux5 <- f[, , it - 1]
    for (l in 1:j) {
      betaprop <- rep(0, p)
      for (t in 1:p) {
        betaprop[t] <- stats::rnorm(1, aux4[l, t], sdpropbeta) * ifelse(l > t, 1, 0) +
          truncnorm::rtruncnorm(1, 0, Inf, aux4[l, t], sdpropbeta2) * ifelse(l == t, 1, 0)
      }

      A1 <- logcondcompbetajProbit(f[, , it - 1], betaprop, alfa[l, ], sigma2[l, it - 1], y, l)

      A2 <- logcondcompbetajProbit(f[, , it - 1], aux4[l, ], alfa[l, ], sigma2[l, it - 1], y, l)
      A <- exp(A1 - A2)
      pa <- min(1, A)

      u <- stats::runif(1)
      if (u < pa) {
        beta[l, , it] <- betaprop
        aux4[l, ] <- betaprop
      } else {
        beta[l, , it] <- beta[l, , it - 1]
      }
    }
    for (l in 1:j) {
      sigma2prop <- truncnorm::rtruncnorm(1, 0, Inf, sigma2[l, it - 1], sdpropsigma2)
      A3 <- logcondcompsigma2jProbit(f[, , it - 1], beta[l, , it], alfa[l, ], sigma2prop, y, l)
      A4 <- logcondcompsigma2jProbit(f[, , it - 1], beta[l, , it], alfa[l, ], sigma2[l, it - 1], y, l)
      B <- exp(A3 - A4)
      pa <- min(1, B)
      u <- stats::runif(1)
      if (u < pa) {
        sigma2[l, it] <- sigma2prop
      } else {
        sigma2[l, it] <- sigma2[l, it - 1]
      }
    }
    for (i in 1:n) {
      fprop <- stats::rnorm(p, aux5[, i], sdpropf)
      A5 <- logcondcompfiProbit(fprop, beta[, , it], alfa, sigma2[, it], y, i)
      A6 <- logcondcompfiProbit(aux5[, i], beta[, , it], alfa, sigma2[, it], y, i)
      C <- exp(A5 - A6)
      pa <- min(1, C)
      u <- stats::runif(1)
      if (u < pa) {
        f[, i, it] <- fprop
        aux5[, i] <- fprop
      } else {
        f[, i, it] <- f[, i, it - 1]
      }
    }
    utils::setTxtProgressBar(pb, it)
  }
  close(pb)
  return(
    list("beta" = beta, "sigma2" = sigma2, "f" = f)
  )
}

fatcat_p_igual_a_1 <- function(y, nit, beta, sigma2, f, alfa, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, n, p, j) {
  pb <- utils::txtProgressBar(min = 2, max = nit, style = 3)
  for (it in 2:nit) {
    aux4 <- as.matrix(beta[, , it - 1])
    aux5 <- t(as.matrix(f[, , it - 1]))
    for (l in 1:j) {
      betaprop <- rep(0, p)
      for (t in 1:p) {
        betaprop[t] <- stats::rnorm(1, aux4[l, t], sdpropbeta) * ifelse(l > t, 1, 0) +
          truncnorm::rtruncnorm(1, 0, Inf, aux4[l, t], sdpropbeta2) * ifelse(l == t, 1, 0)
      }
      A1 <- logcondcompbetajProbit(t(f[, , it - 1]), betaprop, alfa[l, ], sigma2[l, it - 1], y, l)
      A2 <- logcondcompbetajProbit(t(f[, , it - 1]), aux4[l, ], alfa[l, ], sigma2[l, it - 1], y, l)
      A <- exp(A1 - A2)
      pa <- min(1, A)
      u <- stats::runif(1)
      if (u < pa) {
        beta[l, , it] <- betaprop
        aux4[l, ] <- betaprop
      } else {
        beta[l, , it] <- beta[l, , it - 1]
      }
    }
    for (l in 1:j) {
      sigma2prop <- truncnorm::rtruncnorm(1, 0, Inf, sigma2[l, it - 1], sdpropsigma2)
      A3 <- logcondcompsigma2jProbit(t(f[, , it - 1]), as.matrix(beta[l, , it]), alfa[l, ], sigma2prop, y, l)
      A4 <- logcondcompsigma2jProbit(t(f[, , it - 1]), as.matrix(beta[l, , it]), alfa[l, ], sigma2[l, it - 1], y, l)
      B <- exp(A3 - A4)
      pa <- min(1, B)
      u <- stats::runif(1)
      if (u < pa) {
        sigma2[l, it] <- sigma2prop
      } else {
        sigma2[l, it] <- sigma2[l, it - 1]
      }
    }
    for (i in 1:n) {
      fprop <- as.matrix(stats::rnorm(p, aux5[, i], sdpropf))
      A5 <- logcondcompfiProbit(fprop, as.matrix(beta[, , it]), alfa, sigma2[, it], y, i)
      A6 <- logcondcompfiProbit(aux5[, i], as.matrix(beta[, , it]), alfa, sigma2[, it], y, i)
      C <- exp(A5 - A6)
      pa <- min(1, C)
      u <- stats::runif(1)
      if (u < pa) {
        f[, i, it] <- fprop
        aux5[, i] <- fprop
      } else {
        f[, i, it] <- f[, i, it - 1]
      }
    }
    utils::setTxtProgressBar(pb, it)
  }
  close(pb)
  return(
    list("beta" = beta, "sigma2" = sigma2, "f" = f)
  )
}
