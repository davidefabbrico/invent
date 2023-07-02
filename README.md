
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invent

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/treatppmx)](https://CRAN.R-project.org/package=treatppmx)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

An R package for Nonlinear Interactions Varying Coefficient Regression
Model (INteractions Varying coefficiENTs)

The implementation has been done in C++ through the use of Rcpp and
RcppArmadillo.

Authors: Davide Fabbrico and Matteo Pedone

Maintainer: Davide Fabbrico.

## Installation

You can install the development version of invent from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mattpedone/invent")
```

## Usage

In this section, we will demonstrate a basic example to show how the
functions in the R package `invent` work.

The R package contains various functions that are designed to perform
specific tasks. To showcase the functionality, we will go through a
simple example that illustrates the use of one of these functions. This
is what `invent::gendata()` and `invent::invMCMC()` do.

``` r
# Load the R package
library(invent)

# Generate synthetic data
data_list <- gendata(n_obs = 50, p = 10, scenario = 1, nnnc = 3)

# Run the regression model
out <- invMCMC(data_list$Y, data_list$X, iter = 100, burnin = 50, thin = 1)
```
