#' Função para simular dados do modelo categórico
#'
#' @param beta matriz jxp
#' @param Sigma matriz de variancias jxj
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
#' p.real <- 2
#' beta.real <- matrix(
#'   c(
#'     0.99, 0,
#'     0, 0.99,
#'     0.90, 0,
#'     0, 0.90,
#'     0.5, 0.5
#'   ),
#'   nrow = j, ncol = p.real,
#'   byrow = TRUE
#' )
#' Sigma <- diag(c(0.01, 0.05, 0.10, 0.15, 0.20), nrow = j, ncol = j)
#' alfa.real <- matrix(
#'   c(-Inf, qnorm(0.4), qnorm(0.75), qnorm(0.9), Inf),
#'   nrow = j, ncol = (k + 1),
#'   byrow = TRUE
#' )
#' f.real <- t(MASS::mvrnorm(n, rep(0, p.real), diag(1, p.real, p.real)))
#'
#' y <- data_sim(beta.real, Sigma, alfa.real, f.real, link = "probit")
#'
#' @export
#'
data_sim <- function(beta, Sigma, alfa, f, link = c("probit", "logit")) {
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
