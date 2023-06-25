
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fatcat

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/vitorcapdeville/fatcat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vitorcapdeville/fatcat/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The fatcat package provides an interface to fit factor models for
categorical data, described at [Capdeville, Gonçalves and Pereira
(2021)](https://doi.org/10.1080/02664763.2020.1796935).

## Installation

You can install the development version of fatcat from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vitorcapdeville/fatcat")
```

## Usage

``` r
library(fatcat)
#> Loading required package: RcppTN
post = fitfatcat(ordered_data$y, q = 2, nit = 10000,lag = 9, burnin = 1000, quiet = TRUE)
beta_post = posterior::subset_draws(post, variable = "beta", regex = TRUE)
posterior::summarise_draws(beta_post)
#> # A tibble: 10 × 10
#>    variable    mean  median     sd    mad      q5   q95  rhat ess_bulk ess_tail
#>    <chr>      <num>   <num>  <num>  <num>   <num> <num> <num>    <num>    <num>
#>  1 beta_11  1.02    1.02    0.0651 0.0656  0.916  1.13   1.03     50.4     82.7
#>  2 beta_21  0.00780 0.00391 0.0653 0.0638 -0.0970 0.123  1.03     31.5    130. 
#>  3 beta_31  0.966   0.961   0.0684 0.0660  0.857  1.08   1.05     56.6    108. 
#>  4 beta_41  0.0766  0.0723  0.0676 0.0724 -0.0244 0.189  1.06     24.2     94.4
#>  5 beta_51  0.619   0.618   0.0702 0.0703  0.508  0.729  1.04     46.6    162. 
#>  6 beta_12  0       0       0      0       0      0     NA        NA       NA  
#>  7 beta_22  0.936   0.934   0.0651 0.0664  0.833  1.05   1.01    116.     429. 
#>  8 beta_32  0.0973  0.0963  0.0379 0.0366  0.0373 0.166  1.00    229.     380. 
#>  9 beta_42  0.969   0.969   0.0576 0.0545  0.877  1.07   1.02     91.0    179. 
#> 10 beta_52  0.641   0.640   0.0567 0.0576  0.551  0.741  1.01    241.     464.
```

## Documentation

The documentation for the package is available
[here](https://vitorcapdeville.github.io/fatcat/).
