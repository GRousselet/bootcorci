% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{twocorci}
\alias{twocorci}
\title{Compare two independent correlations}
\usage{
twocorci(
  x1,
  y1,
  x2,
  y2,
  method = "pbcor",
  nboot = 2000,
  alpha = 0.05,
  alternative = "two.sided",
  null.value = 0,
  saveboot = TRUE,
  ...
)
}
\arguments{
\item{x1, y1}{Two dependent vectors of the same length from group 1.}

\item{x2, y2}{Two dependent vectors of the same length from group 2.}

\item{method}{A function that returns a correlation coefficient.
Options in \code{bootcorci} include "pearson", "spearman", "pbcor", "wincor".
Default is "pbcor".}

\item{nboot}{Number of bootstrap samples. Default 2000.}

\item{alpha}{Alpha level. Default 0.05. For corfun = pearson, alpha is restricted
to 0.05 because the confidence interval adjustments have not been calculated for
other alphas.}

\item{alternative}{Type of test, either "two.sided" (default), "greater" for
positive correlations, or "less" for negative correlations.}

\item{null.value}{Hypothesis to test. Default 0.}

\item{saveboot}{Option to save bootstrap samples. Default TRUE.}

\item{...}{Optional parameter to pass to correlation function.}
}
\value{
\itemize{
\item \code{estimate1} the correlation coefficient for group 1.
\item \code{estimate2} the correlation coefficient for group 2.
\item \code{difference} the difference between correlation coefficients.
\item \code{conf.int} the bootstrap confidence interval.
\item \code{p.value} the bootstrap p value.
\item \code{bootsamples} the bootstrap samples.
}
}
\description{
Compute a 1-alpha percentile bootstrap confidence interval for the difference between
two correlation coefficients corresponding to two independent groups.
The default correlation is the percentage bend.
When using Pearson's correlation, the confidence interval is adjusted
to compensate for the error term's heteroscedasticity.
Missing values are automatically removed.
}
\section{Note}{
 Modified from functions \code{twocor} and \code{twopcor}
from Rallfun-v37.txt - see \url{https://github.com/nicebread/WRS/} and
\url{http://dornsife.usc.edu/labs/rwilcox/software/}.
}

\references{
Wilcox, R.R. (2009) Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality.
Communications in Statistics - Simulation and Computation, 38, 2220–2234.

Wilcox, R.R. (2017) Introduction to Robust Estimation and Hypothesis Testing, 4th edition. Academic Press.
}
