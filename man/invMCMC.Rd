% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SamplerCpp.R
\name{invMCMC}
\alias{invMCMC}
\title{invMCMC}
\usage{
invMCMC(
  y,
  x,
  y_val = NULL,
  x_val = NULL,
  hyperpar = c(3, 1, 1, 1, 0.00025, 0.4, 1.6, 0.2, 1.8, 0.4, 1.6, 0.2, 1.8),
  mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5),
  rank = 0.95,
  iter = 10000,
  burnin = iter/2,
  thin = 5,
  ha = 2,
  detail = FALSE,
  data = NULL,
  pb = TRUE
)
}
\description{
invMCMC
}
