
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invent

<!-- badges: start -->

[![CRAN
status](https://img.shields.io/cran/v/invent)](https://CRAN.R-project.org/package=invent)
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
devtools::install_github("davidefabbrico/invent")
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
data_list <- invent::gendata(n_obs = 400, p = 10, scenario = 4)

y <- data_list$Y
x <- data_list$X

# Run the regression model
out <- invent::invMCMC(y, x, iter = 10000, burnin = 5000, thin = 5, ha = 2, 
                       data = data_list, detail = TRUE, pb = TRUE)
```

For additional details, please refer to the howtouse (.Rmd or .html) file available in this repository.
