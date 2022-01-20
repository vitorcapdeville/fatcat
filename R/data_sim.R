#' Função para simular dados do modelo categórico
#'
#' @param beta matriz jxp
#' @param sigma vetor com a diagonal da matriz de variancias, tamanho j
#' @param alfa matriz j x (k+1)
#' @param f matriz pxn
#' @param link probit ou logit
#'
#' @return matriz jxn
#'
#' @examples
#' j <- 5
#' n <- 300
#' k <- 4
#' true_p <- 2
#'
#' true_beta <- matrix(
#'   c(
#'     0.99, 0,
#'     0, 0.99,
#'     0.90, 0,
#'     0, 0.90,
#'     0.5, 0.5
#'   ),
#'   nrow = j, ncol = true_p,
#'   byrow = TRUE
#' )
#'
#' true_sigma <- c(0.01, 0.05, 0.10, 0.15, 0.20)
#'
#' true_alfa <- matrix(
#'   c(-Inf, qnorm(0.4), qnorm(0.75), qnorm(0.9), Inf),
#'   nrow = j, ncol = (k + 1),
#'   byrow = TRUE
#' )
#'
#' true_f <- t(
#'   MASS::mvrnorm(n, rep(0, true_p), diag(1, true_p, true_p))
#' )
#'
#' y <- data_sim(true_beta, true_sigma, true_alfa, true_f, link = "probit")
#'
#' @export
#'
data_sim <- function(beta, sigma, alfa, f, link = c("probit", "logit")) {
  Sigma <- diag(sigma)
  stopifnot(nrow(beta) == nrow(alfa))
  stopifnot(ncol(beta) == nrow(f))
  stopifnot(nrow(Sigma) == ncol(Sigma))
  stopifnot(nrow(beta) == ncol(Sigma))
  j <- nrow(beta)
  n <- ncol(f)
  k <- ncol(alfa) - 1

  # Para a geracao logit
  s <- sqrt(diag(Sigma) * 3 / (pi^2))

  y_star <- matrix(NA, nrow = j, ncol = n)
  E <- matrix(NA, nrow = j, ncol = n)
  for (i in 1:n) {
    E[, i] <- MASS::mvrnorm(1, rep(0, j), Sigma)
    for (l in 1:j) {
      y_star[l, i] <- beta[l, ] %*% as.matrix(f[, i]) + ifelse(link == "probit", E[l, i], stats::rlogis(1, 0, s[l]))
    }
  }

  y <- array(0, dim = c(j, n))
  for (l in 1:j) {
    for (i in 1:n) {
      for (s in 1:k) {
        if (y_star[l, i] >= alfa[l, s] && y_star[l, i] < alfa[l, s + 1]) {
          y[l, i] <- s
        }
      }
    }
  }
  return(y)
}
