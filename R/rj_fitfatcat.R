# Criar funcao que usa ou o rjmcmc ou o fitfatcat
# if (!is.null(p)) {
#   message("Usando p como fixo e igual a ", p, ". Para estimar p via RJMCMC, coloque p = NULL.")
#   return(fitfatcat(y = y, p = p, nit = nit, ...))
# } else {


#' Funcao auxiliar para identificar o maior valor possivel do numero de fatores
#' para um dado numero de variáveis.
#'
#' @param m numero de variaveis, nrow(y)
#'
#' @return limete superior para o numero de fatores
#'
#' @export
q_max <- function(m) {
  k <- seq(m)
  aux <- (m * (m + 1) / 2 - m * (k + 1) + k * (k - 1) / 2)
  return(max(which(aux >= 0)))
}

#' Calculo da moda
#'
#' @param x vetor
#'
#' @return moda do vetor
#'
Mode <- function(x) {
  d <- stats::density(x)
  aux <- which.max(d$y)
  return(d$x[aux])
}

#' Algoritmo de RJMCMC para o modelo fatorial categórico para dados ordenados
#'
#' @param y dados. deve ter dimensao jxn
#' @param alpha pontos de corte. usualmente fixados no valor verdadeiro
#' ou na estimativa usando a correlacao policorica.
#' @param chains cadeias com amostras da posteriori dos parametros dado p,
#' para todos os p's possiveis.
#' @param nit_rjmcmc numero de iteracoes para o rjmcmc
#' @param init_q numero inicial para o numero de fatores
#' @param C0,a,b hiperparametros das prioris do MCMC. Devem coincidir com as cadeias geradas em chains.
#' @param a_rj,b_rj hiper parametros para a proposta do RJMCMC
#' @param sdpropbeta,sdpropbeta2,sdpropf,sdpropsigma2 desvio padrao das propostas.
#' @param ... argumentos extras.
#'
#' @export
#'
rjmcmc_fitfatcat <- function(y, alpha, chains, nit_rjmcmc = 500, init_q = 2, C0 = 10000, a = 0.001, b = 0.001, a_rj = 2, b_rj = 1, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075, ...) {
  p <- nrow(y)
  n <- ncol(y)

  nit <- dim(chains[[1]]$beta)[3]

  max_q <- q_max(p)
  stopifnot(max_q == length(chains))

  # Calcule parametros da proposta pro rjmcmc
  b_k <- purrr::map(
    chains,
    function(x) {
      matrix(apply(x$beta, seq(length(dim(x$beta)) - 1), mean), nrow = dim(x$beta)[1])
    }
  )

  B_k <- purrr::map(
    chains,
    function(x) {
      lapply(seq(dim(x$f)[1]), function(y) stats::cov(t(x$beta[, y, ])))
    }
  )

  v <- purrr::map(chains, function(x) apply(x$sigma2, 1, Mode))

  # Definir onde as cadeias do RJMCMC ficarao armazenadas
  beta_final <- array(NA, dim = c(p, max_q, nit_rjmcmc))
  sigma2_final <- array(NA, dim = c(p, nit_rjmcmc))
  f_final <- array(NA, dim = c(max_q, n, nit_rjmcmc))
  q_final <- rep(NA, nit_rjmcmc)

  # Valores iniciais (o final da cadeia ja aquecida.) [4.2, 0]
  beta <- array(0, dim = c(p, init_q))
  sigma2 <- array(NA, dim = c(p))
  f <- array(NA, dim = c(init_q, n))
  beta <- chains[[init_q]]$beta[, , nit]
  sigma2 <- chains[[init_q]]$sigma2[, nit]
  f <- chains[[init_q]]$f[, , nit]

  curr_theta <- algoritmo_probit(
    y, 2, beta, f, sigma2, alpha,
    C0, a, b,
    sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2,
    F
  )

  curr_beta <- matrix(curr_theta$beta, nrow = nrow(curr_theta$beta))
  curr_sigma2 <- curr_theta$sigma2
  curr_f <- matrix(curr_theta$beta, nrow = nrow(curr_theta$beta))

  beta_final[, 1:init_q, 1] <- matrix(curr_theta$beta, nrow = nrow(curr_theta$beta))
  sigma2_final[, 1] <- curr_theta$sigma2
  f_final[1:init_q, , 1] <- matrix(curr_theta$f, nrow = nrow(curr_theta$f))
  q_final[1] <- init_q

  pb <- utils::txtProgressBar(min = 1, max = nit_rjmcmc, initial = 1, style = 3)
  for (it in 2:nit_rjmcmc) {
    # Defina 1 como subir 1 dimensao
    # Defina -1 como descer uma dimensao
    # Defina a probabilidade de subir de dimensao como 1/2

    curr_q <- q_final[it - 1]
    curr_beta <- matrix(beta_final[, 1:curr_q, it - 1], nrow = p)
    curr_sigma2 <- sigma2_final[, it - 1]
    curr_f <- matrix(f_final[1:curr_q, , it - 1], ncol = n)

    if (curr_q == 1) {
      mov <- 1
    } else if (curr_q == max_q) {
      mov <- -1
    } else {
      mov_coin <- stats::runif(1)
      if (mov_coin > 0.5) {
        mov <- 1
      } else {
        mov <- -1
      }
    }

    # proposta de mudanca de dimensao
    prop_q <- curr_q + mov

    # amostra de beta dada a nova dimensao
    prop_beta <- matrix(0, nrow = p, ncol = prop_q)
    for (l in 1:prop_q) {
      prop_beta[, l] <- MASS::mvrnorm(1, b_k[[prop_q]][, l], b_rj * B_k[[prop_q]][[l]])
    }

    # amostra de sigma2 dada a nova dimensao
    prop_sigma2 <- rep(0, p)
    for (j in 1:p) {
      prop_sigma2[j] <- invgamma::rinvgamma(1, shape = a_rj, rate = a_rj * v[[prop_q]][j])
    }

    # 3 - Calcular a taxa de aceitacao/rejeicao e escolhe aceitar ou rejeitar
    # o pulo para o k proposto.

    log_num <- 0
    log_den <- 0

    # Verossimilhanca e priori de beta, no valor proposto
    log_num <- l_vero_marginal_probit(y, alpha, prop_beta, prop_sigma2) +
      l_priori_beta(prop_beta, C0) + l_priori_sigma2(prop_sigma2, a, b)
    # Priori uniforme para o numero de fatores, logo, nao precisa de nada aqui
    # Probabilidades de transicao p -> p' e p' -> p sao sempre iguais.

    # Densidade da proposta no valor anterior
    log_d_prop_beta <- 0
    for (l in 1:curr_q) {
      zeros <- which(curr_beta[, l] == 0)
      if (length(zeros) > 0) {
        aux_curr_beta <- curr_beta[, l][-zeros]
        aux_b_k <- b_k[[curr_q]][, l][-zeros]
        aux_B_k <- B_k[[curr_q]][[l]][-zeros, -zeros]
      } else {
        aux_curr_beta <- curr_beta[, l]
        aux_b_k <- b_k[[curr_q]][, l]
        aux_B_k <- B_k[[curr_q]][[l]]
      }
      log_d_prop_beta <- log_d_prop_beta + mvtnorm::dmvnorm(aux_curr_beta, aux_b_k, b_rj * aux_B_k, log = T)
    }

    log_num <- log_num + sum(invgamma::dinvgamma(curr_sigma2, shape = a_rj, rate = a_rj * v[[curr_q]], log = T)) + log_d_prop_beta

    # Verossimilhanca e priori de beta, no valor anterior
    log_den <- l_vero_marginal_probit(y, alpha, curr_beta, curr_sigma2) +
      l_priori_beta(curr_beta, C0) + l_priori_sigma2(curr_sigma2, a, b)

    # Densidade da proposta no valor proposto
    log_d_prop_beta <- 0
    for (l in 1:prop_q) {
      zeros <- which(prop_beta[, l] == 0)
      if (length(zeros) > 0) {
        aux_prop_beta <- prop_beta[, l][-zeros]
        aux_b_k <- b_k[[prop_q]][, l][-zeros]
        aux_B_k <- B_k[[prop_q]][[l]][-zeros, -zeros]
      } else {
        aux_prop_beta <- prop_beta[, l]
        aux_b_k <- b_k[[prop_q]][, l]
        aux_B_k <- B_k[[prop_q]][[l]]
      }
      log_d_prop_beta <- log_d_prop_beta + mvtnorm::dmvnorm(aux_prop_beta, aux_b_k, b_rj * aux_B_k, log = T)
    }
    log_den <- log_den + sum(invgamma::dinvgamma(prop_sigma2, shape = a_rj, rate = a_rj * v[[prop_q]], log = T)) + log_d_prop_beta

    log_ratio <- log_num - log_den
    rho <- exp(log_ratio)
    coin_accept <- stats::runif(1)
    alpha_mov <- min(1, rho)
    if (is.na(alpha_mov)) browser()
    if (coin_accept < alpha_mov) {
      # Se k é aceito, salva k e theta k.
      # Aceitamos a proposta de mover dimensao
      q_final[it] <- prop_q
      # O que entra no valor de f aqui? Ele também depende de k...
      # Aqui deveria ser marginal em f, talvez? Mas ai como eu geraria uma amostra de f?
      # Movemos para a nova dimensao
      curr_theta <- algoritmo_probit(
        y, 2, matrix(prop_beta, nrow = p), matrix(f_final[1:prop_q, , it - 1], ncol = n), prop_sigma2, alpha,
        C0, a, b,
        sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2,
        F
      )

      beta_final[, 1:prop_q, it] <- curr_theta$beta
      sigma2_final[, it] <- curr_theta$sigma2
      f_final[1:prop_q, , it] <- curr_theta$f
    } else {
      # Se k é rejeitado, usa o k anterior e theta k anterior.
      q_final[it] <- curr_q
      # Rejeitamos a mudanca de dimensao
      # Amostramos um valor da dimensao corrente.

      curr_theta <- algoritmo_probit(
        y, 2, curr_beta, curr_f, curr_sigma2, alpha,
        C0, a, b,
        sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2,
        F
      )

      beta_final[, 1:curr_q, it] <- curr_theta$beta
      sigma2_final[, it] <- curr_theta$sigma2
      f_final[1:curr_q, , it] <- curr_theta$f
    }
    utils::setTxtProgressBar(pb, it)
  }
  return(list(beta = beta_final, sigma2 = sigma2_final, f = f_final, q = q_final))
}
