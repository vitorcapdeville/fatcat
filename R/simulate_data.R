#' Simulate data from a factor model for categorical ordinal data.
#'
#'
#' The simulated data follows the model:
#'
#' y_i|f_i = beta f_i + epsilon_i
#' f_i ~ N(0, I_q)
#'
#' where:
#' y_i is a vector of dimension j
#' beta is a matrix of dimension j x q
#' f_i is a vector of dimension q
#' epsilon_i is a vector of dimension j
#' q is the number of factors
#' j is the number of variables
#'
#' The output of the simulation is y_i^*, which is a categorization of y_i.
#' y_ij^* = k if alpha_{k-1} <= y_ij <= alpha_k
#' The number of categories is fixed at k. There are k + 1 cutoff points,
#' and usually the first cutoff point is -infty and the last one is +infty
#'
#' The model can also be written as:
#'
#' y = beta beta' + Sigma
#'
#' which is known as the marginal model.
#'
#' The second form of the model can be simulated from by providing `n` instead of `f`.
#' For now, simulating from the marignal model is compatible only with the probit link function.
#'
#'
#' @param beta A matrix with dimension j x q, where j is the number of variables and q the number of latent factors.
#' @param sigma A vector with the elements of the diagonal of the matrix of variances, with j elements.
#' @param alpha A matrix of dimension j x (k+1).
#' @param f A matrix of dimension p x n. Should only be provided if `n` is not provided.
#' @param n The number of observations to be simulated. Should only be provided if `f` is not provided.
#' @param link The link function to be used, either "probit" or "logit".
#'
#' @return A matrix with dimension j x n, with the simulated data.
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
#' true_alpha <- matrix(
#'   c(-Inf, qnorm(0.4), qnorm(0.75), qnorm(0.9), Inf),
#'   nrow = j, ncol = (k + 1),
#'   byrow = TRUE
#' )
#'
#' true_f <- t(
#'   MASS::mvrnorm(n, rep(0, true_p), diag(1, true_p, true_p))
#' )
#' y <- simulate_data(
#'   beta = true_beta,
#'   sigma = true_sigma,
#'   alpha = true_alpha,
#'   f = true_f,
#'   link = "probit"
#' )
#'
#' y_marginal <- simulate_data(
#'   beta = true_beta,
#'   sigma = true_sigma,
#'   alpha = true_alpha,
#'   n = n,
#'   link = "probit"
#' )
#'
#' @export
#'
simulate_data <- function(beta, sigma, alpha, n = NULL, f = NULL, link = c("probit", "logit")) {
  link = match.arg(link)
  Sigma <- diag(sigma)
  marginal = is.null(f)

  if (is.null(f) & is.null(n)) stop("Either `n` or `f` must be provided.")
  if (!is.null(f) & !is.null(n)) stop("One and only one of `n` and `f` must be provided.")

  stopifnot(nrow(beta) == nrow(alpha))
  stopifnot(ncol(beta) == nrow(f))
  stopifnot(nrow(Sigma) == ncol(Sigma))
  stopifnot(nrow(beta) == ncol(Sigma))
  j <- nrow(beta)
  if (!marginal) n <- ncol(f)
  k <- ncol(alpha) - 1

  # Para a geracao logit
  s <- sqrt(diag(Sigma) * 3 / (pi^2))

  y_star <- matrix(NA, nrow = j, ncol = n)
  E <- matrix(NA, nrow = j, ncol = n)
  if (marginal) {
    stopifnot(link == "probit") # Marginal feito apenas para probit por enquanto
    for (i in 1:n) {
      y_star[, i] <- MASS::mvrnorm(1, rep(0, j), beta %*% t(beta) + Sigma)
    }
  }else{
    for (i in 1:n) {
      E[, i] <- MASS::mvrnorm(1, rep(0, j), Sigma)
      for (l in 1:j) {
        y_star[l, i] <- beta[l, ] %*% as.matrix(f[, i]) + {if (link == "probit") E[l, i] else stats::rlogis(1, 0, s[l])}
      }
    }
  }

  y <- array(0, dim = c(j, n))
  for (l in 1:j) {
    for (i in 1:n) {
      for (s in 1:k) {
        if (y_star[l, i] >= alpha[l, s] && y_star[l, i] < alpha[l, s + 1]) {
          y[l, i] <- s
        }
      }
    }
  }
  return(y)
}
