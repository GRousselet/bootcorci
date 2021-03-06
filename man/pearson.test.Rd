% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation.R
\name{pearson.test}
\alias{pearson.test}
\title{Compute Pearson's rho}
\usage{
pearson.test(x, y, alternative = "two.sided")
}
\arguments{
\item{x, y}{Two vectors of the same length.}
}
\value{
\itemize{
  \item \code{estimate} the correlation coefficient.
  \item \code{p.value} the p-value of the test.
  \item \code{statistic} the t statistic of the test.
  }
}
\description{
Missing values are automatically removed.
}
