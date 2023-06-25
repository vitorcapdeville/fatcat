array2matrix <- function(array, name) {
  new_as_matrix <- matrix(array, nrow = dim(array)[3], byrow = T)
  aux_names <- expand.grid(seq_len(dim(array)[1]), seq_len(dim(array)[2]))
  colnames(new_as_matrix) <- paste0(name, "_", aux_names$Var1, aux_names$Var2)
  return(new_as_matrix)
}


#' Draws samples from the posterior of the factorial model for categorical data.
#'
#' @param y A matrix of integers with dimension j x n, where j is the number of variables and n the number of observations.
#' @param q An integer with the number of factors.
#' @param nit An integer with the number of iterations.
#' @param burnin An integer with the number of initial iterations to be ignored.
#' @param lag An integer with the lag of the iterations taken. ``lag = 1`` means take every iteration. ``lag = n`` means
#'  take one iteration every n iterations.
#' @param init_beta An integer with the initial value of the first column of beta.
#' @param init_sigma2 An integer with the initial value of sigma2. All elements of the sigma2 vector will be initialized
#'  with the same value.
#' @param C0 An integer with the variance of the prior for beta.
#' @param a,b An integer with the shape and scale parameters of the prior of sigma2.
#' @param sdpropbeta An integer with the standard deviation of the proposal of betas that are unrestricted.
#' @param sdpropbeta2 An integer with the standard deviation of the proposal of betas that are restricted to be greater
#'  than or equal to zero.
#' @param sdpropf An integer with the standard deviation of the proposal of factors.
#' @param sdpropsigma2  An integer with the standard deviation of the proposal of sigma2.
#' @param dist An string with the name of the error distribution. Currently, only `probit` is accepted.
#' @param quiet If TRUE, prevents the display of progress bar.
#' @param alpha An vector with the true value of alpha. Usually, starts with -Inf and ends with +Inf.
#'
#' @return A posterior object with a sample from the posterior distribution.
#'
#' @examples
#' post = fitfatcat(ordered_data$y, q = 2, nit = 1000,lag = 9, burnin = 100, quiet = TRUE)
#' beta_post = posterior::subset_draws(post, variable = "beta", regex = TRUE)
#' posterior::summarise_draws(beta_post)
#'
#' @import RcppTN
#'
#' @export
fitfatcat <- function(y, q, nit, burnin = 0, lag = 1, init_beta = 0.5, init_sigma2 = 2, C0 = 100000, a = 0.001, b = 0.001, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.075, dist = c("probit"), quiet = F, alpha = NULL) {
  # Definicoes e valores iniciais
  dist <- match.arg(dist)

  p <- nrow(y)
  n <- ncol(y)

  beta <- array(0, dim = c(p, q))
  sigma2 <- array(NA, dim = c(p))
  f <- array(NA, dim = c(q, n))
  beta[, 1] <- init_beta
  sigma2 <- rep(init_sigma2, p)
  f <- matrix(0.5, nrow = q, ncol = n)

  if (is.null(alpha)) alpha <- cbind(-Inf, psych::polychoric(t(y))$tau, Inf)

  # Algoritmo
  res <- algoritmo_probit(y, nit + 1, beta, f, sigma2, alpha, C0, a, b, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, !quiet)

  new_beta <- array2matrix(res$beta, "beta")
  new_f <- array2matrix(res$f, "f")
  new_sigma2 <- t(res$sigma2)
  colnames(new_sigma2) <- paste0("sigma2_", seq_len(dim(res$sigma2)[1]))
  posterior_dist <- posterior::as_draws(cbind(new_beta, new_f, new_sigma2))

  posterior_dist <- posterior::subset_draws(posterior_dist, iteration = seq(from = burnin + 1, to = min(nit, dim(res$beta)[3])))
  posterior_dist <- posterior::thin_draws(posterior_dist, thin = lag)


  return(posterior_dist)
}
