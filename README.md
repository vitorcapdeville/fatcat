
# fatcat

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The fatcat package provides an interface to fit factor models for categorical data, described at [Capdeville, Gonçalves and Pereira (2021)]( https://doi.org/10.1080/02664763.2020.1796935).

## Installation

You can install the development version of fatcat from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vitorcapdeville/fatcat")
```

## Usage

``` r
library(fatcat)
post = fitfatcat(ordered_data$y, q = 2, nit = 10000,lag = 9, burnin = 1000)
beta_post = posterior::subset_draws(res, variable = "beta", regex=T)
posterior::summarise_draws(beta_post)
```

## Documentation

The documentation for the package is available [here](https://vitorcapdeville.github.io/fatcat/).

