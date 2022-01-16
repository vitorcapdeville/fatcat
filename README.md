
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fatcat

<!-- badges: start -->
<!-- badges: end -->

Fit bayesian factor models for categorical data.

## Installation

You can install the development version of fatcat like so:

``` r
install.packages(remotes)
remotes::install_github("vitorcapdeville/fatcat")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fatcat)
# Genarate some fake data
j <- 5
n <- 300
k <- 4
p.real <- 2
beta.real <- matrix(
  c(
    0.99, 0,
    0, 0.99,
    0.90, 0,
    0, 0.90,
    0.5, 0.5
  ),
  nrow = j, ncol = p.real,
  byrow = TRUE
)
Sigma <- diag(c(0.01, 0.05, 0.10, 0.15, 0.20), nrow = j, ncol = j)
alfa.real <- matrix(c(-Inf, qnorm(0.4), qnorm(0.75), qnorm(0.9), Inf), nrow = j, ncol = (k + 1), byrow = T)
f.real <- t(MASS::mvrnorm(n, rep(0, p.real), diag(1, p.real, p.real)))

y <- data_sim(beta.real, Sigma, alfa.real, f.real, link = "probit")

# Visualizing the variance structure.
cor <- psych::polychoric(t(y))$rho
corrplot::corrplot(cor, type = "upper")
# Note that variables with a high beta for the same factor have high correlation
# with each other.

# Fit the model
nit <- 1000
res <- fitfatcat(y, p = 2, nit)

# See the results
require(ggplot2)
require(patchwork)

plot_beta <- function(indice, beta, beta.real) {
  ggplot(mapping = aes(x = 1:dim(beta)[3])) +
    ylab(glue::glue("beta[", indice[1], ",", indice[2], "]")) +
    xlab("") +
    theme_bw() +
    geom_line(aes(y = res$beta[indice[1], indice[2], ])) +
    geom_line(aes(y = beta.real[indice[1], indice[2]]), color = "red")
}

indices <- matrix(unlist(expand.grid(1:j, 1:p.real)), nrow = j*2, ncol = p.real, byrow = F)

plots_beta <- apply(indices, 1, plot_beta, beta = res$beta, beta.real = beta.real)

wrap_plots(plots_beta, nrow = j, ncol = p.real)
```
