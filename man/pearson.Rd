% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_simplified.R
\name{pearson}
\alias{pearson}
\title{Compute Pearson's rho}
\usage{
pearson(x, y)
}
\arguments{
\item{x, y}{Two vectors of the same length.}
}
\value{
The correlation coefficient.
}
\description{
Simplified version of \code{pearson.test} which returns only the correlation coefficient.
Called by the bootstrap functions for faster execution.
}
