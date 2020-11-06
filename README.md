
<!-- badges: start -->
[![R build status](https://github.com/mrc-ide/markovid/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/markovid/actions)
<!-- badges: end -->

# markovid

This is an MCMC package for carrying out preliminary fits of hospital progression parameters (transition probabilities and durations in each state). The outputs of these fits can then be passed up to the larger mechanistic model as priors or relative risks for a further level of fitting.

NB, this package is essentially a re-purpose of the drjacoby package, and still contains some legacy code (for example objects of class "drjacoby").

The R_ignore directory contains two scripts that are useful for getting a feel for the package. The first demonstrates how to simulate data from the underlying progression model, and the second demonstrates how to run the MCMC on example data.
