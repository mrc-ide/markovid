% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R
\name{get_spline}
\alias{get_spline}
\title{Get cubic splines from posterior draws}
\usage{
get_spline(mcmc_samples, nodex, age_vec, scale = 1)
}
\arguments{
\item{mcmc_samples}{matrix of posterior draws.}

\item{nodex}{x-coordinates of nodes.}

\item{age_vec}{x-coordinates at which to evaluate cubic splines.}

\item{scale}{posterior splines are logistic transformed and scaled to the
interval \code{[0,scale]}.}
}
\description{
Given a matrix of posterior draws of node locations, calculate
  cubic splines over all draws, then apply logistic transform to get
  transformed spline in a finite range.
}
