## code to prepare `ordered_data` dataset goes here
library(fatcat)
set.seed(1)

j <- 5
n <- 300
k <- 4
true_p <- 2

true_beta <- matrix(
  c(
    0.99, 0.00,
    0.00, 0.99,
    0.90, 0.00,
    0.00, 0.90,
    0.50, 0.50
  ),
  nrow = j, ncol = true_p,
  byrow = TRUE
)

true_sigma <- c(0.01, 0.05, 0.10, 0.15, 0.20)

true_alpha <- matrix(
  c(-Inf, qnorm(0.4), qnorm(0.75), qnorm(0.9), Inf),
  nrow = j, ncol = (k + 1),
  byrow = TRUE
)

y <- simulate_data(beta = true_beta, sigma = true_sigma, alpha = true_alpha, n = n, link = "probit")

ordered_data <- list(y = y, true_beta = true_beta, true_sigma = true_sigma, true_alpha = true_alpha)

usethis::use_data(ordered_data, overwrite = TRUE)
