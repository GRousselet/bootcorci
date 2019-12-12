#' Compute a 1-alpha percentile bootstrap confidence interval for a correlation.
#'
#' The default correlation is the percentage bend correlation.
#' When using Pearson's correlation, the confidence interval is adjusted
#' to compensate for the error term's heteroscedasticity.
#' Missing values are automatically removed.
#'
#' @param x,y Two vectors of the same length.
#' @param method A function that returns a correlation coefficient.
#' Options in \code{bootcorci} include "pearson", "spearman", "pbcor", "wincor".
#' Default is "pbcor".
#' @param nboot Number of bootstrap samples. Default 2000.
#' @param alpha Alpha level. Default 0.05. For corfun = pearson, alpha is restricted
#' to 0.05 because the confidence interval adjustments have not been calculated for
#' other alphas.
#' @param alternative Type of test, either "two.sided", "greater" for positive
#'   correlations, or "less" for negative correlations.
#' @param null.value Hypothesis to test. Default 0.
#' @param ... Optional parameter to pass to correlation function.
#'
#' @return \itemize{
#' \item \code{estimate} the correlation coefficient.
#' \item \code{p.value} the p-value of the test.
#' \item \code{statistic} the t statistic of the test.
#' \item \code{bootsamples} the bootstrap samples.
#' }
#'
#' @section Note: Modified from functions \code{corb} and \code{pcorb}
#' from Rallfun-v37.txt - see \url{https://github.com/nicebread/WRS/} and
#' \url{http://dornsife.usc.edu/labs/rwilcox/software/}.
#'
#' @export
corci <- function(x, y, method="pbcor", nboot=2000, alpha=.05,
                  alternative = "two.sided", null.value=0, ...){
  if(!is.vector(x) || !is.vector(y)){
    stop("corci: x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("corci: the vectors do not have equal lengths.")
  }
  m <- cbind(x,y)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x <- m[,1]
  y <- m[,2]
  eval(parse(text=paste("corfun=",method)))
  est <- corfun(x, y, ...)
  data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
  bvec <- apply(data, 1, corbsub, x, y, corfun, ...) # Create a 1 by nboot matrix.
  if(method == "pearson"){ # Get quantile indices
    out <- adj.corbootci(n, nboot)
    ilow <- out[1]
    ihi <- out[2]
  } else {
    ilow <- floor((alpha/2)*nboot+.5)
    ihi <- floor((1-alpha/2)*nboot+.5)
  }
  bsort <- sort(bvec)
  corci <- c(bsort[ilow], bsort[ihi])
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==null.value))/nboot
    sig <- 2 * min(phat, 1 - phat)
  }
  if(alternative == "greater"){
    sig <- 1 - sum(bvec >= null.value)/nboot
  }
  if(alternative == "less"){
    sig <- 1 - sum(bvec <= null.value)/nboot
  }

  list(conf.int=corci, p.value=sig, estimate=est, bootsamples=bvec)
}

# ::::::::::::::::::::::::::::

#' Compare two independent correlations
#'
#' Compute a 1-alpha percentile bootstrap confidence interval for the difference between
#' two correlation coefficients corresponding to two independent groups.
#' The default correlation is the percentage bend.
#' When using Pearson's correlation, the confidence interval is adjusted
#' to compensate for the error term's heteroscedasticity.
#' Missing values are automatically removed.
#'
#' @param x1,y1 Two dependent vectors of the same length from group 1.
#' @param x2,y2 Two dependent vectors of the same length from group 2.
#' @param method A function that returns a correlation coefficient.
#' Options in \code{bootcorci} include "pearson", "spearman", "pbcor", "wincor".
#' Default is "pbcor".
#' @param nboot Number of bootstrap samples. Default 2000.
#' @param alpha Alpha level. Default 0.05. For corfun = pearson, alpha is restricted
#' to 0.05 because the confidence interval adjustments have not been calculated for
#' other alphas.
#' @param alternative Type of test, either "two.sided" (default), "greater" for
#' positive correlations, or "less" for negative correlations.
#' @param null.value Hypothesis to test. Default 0.
#' @param ... Optional parameter to pass to correlation function.
#'
#' @return \itemize{
#' \item \code{estimate1} the correlation coefficient for group 1.
#' \item \code{estimate2} the correlation coefficient for group 2.
#' \item \code{difference} the difference between correlation coefficients.
#' \item \code{conf.int} the bootstrap confidence interval.
#' \item \code{p.value} the bootstrap p value.
#' \item \code{bootsamples} the bootstrap samples.
#' }
#'
#' @section Note: Modified from functions \code{twocor} and \code{twopcor}
#' from Rallfun-v37.txt - see \url{https://github.com/nicebread/WRS/} and
#' \url{http://dornsife.usc.edu/labs/rwilcox/software/}.
#'
#' @export
twocorci <- function(x1, y1, x2, y2,
                     method="pbcor", nboot=2000, alpha=.05,
                     alternative = "two.sided", null.value=0, ...){
  if(!is.vector(x1) || !is.vector(x2) || !is.vector(y1) || !is.vector(y2)){
    stop("twocorci: x1, x2, y1, y2  must be vectors.")
  }
  if(length(x1)!=length(y1)){
    stop("twocorci: x1 and y1 do not have equal lengths.")
  }
  if(length(x2)!=length(y2)){
    stop("twocorci: x2 and y2 do not have equal lengths.")
  }
  m <- cbind(x1,y1)
  m <- m[complete.cases(m), ]
  n1 <- nrow(m)
  x1 <- m[,1]
  y1 <- m[,2]
  m <- cbind(x2,y2)
  m <- m[complete.cases(m), ]
  n2 <- nrow(m)
  x2 <- m[,1]
  y2 <- m[,2]
  # correlations
  eval(parse(text=paste("corfun=",method)))
  r1 <- corfun(x1,y1,...)
  r2 <- corfun(x2,y2,...)
  # bootstrap
  data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
  bvec1 <- apply(data1, 1, corbsub, x1, y1, corfun, ...) # A 1 by nboot matrix.
  data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
  bvec2 <- apply(data2, 1, corbsub, x2, y2, corfun, ...) # A 1 by nboot matrix.
  bvec <- bvec1 - bvec2
  bsort <- sort(bvec)
  if(method == "pearson"){ # Get quantile indices
    out <- adj.corbootci(n, nboot)
    ilow <- out[1]
    ihi <- out[2]
  } else {
    ilow <- floor((alpha/2)*nboot+.5)
    ihi <- floor((1-alpha/2)*nboot+.5)
  }
  corci <- c(bsort[ilow], bsort[ihi])
  # p value
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==0))/nboot
    sig <- 2 * min(phat, 1 - phat)
  }
  if(alternative == "greater"){
    sig <- 1 - sum(bvec >= null.value)/nboot
  }
  if(alternative == "less"){
    sig <- 1 - sum(bvec <= null.value)/nboot
  }
  list(estimate1=r1, estimate2=r2, difference=r1-r2,
       conf.int=corci, p.value=sig, bootsamples=bvec)
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::

#' Compare two dependent correlations: overlapping case
#'
#' Compute a 1-alpha percentile bootstrap confidence interval for the difference
#' between two correlation coefficients corresponding to two dependent groups,
#' in the overlapping case. Compare correlation between x1 and y to the
#' correlation between x2 and y. The default correlation is the percentage
#' bend correlation. This function is inappropriate to make inferences about
#' Pearson's correlations. Missing values are automatically removed.
#'
#' @param x1,x2,y Three dependent vectors of the same length.
#' @param method A function that returns a correlation coefficient.
#' Options include "spearman", "pbcor", "wincor". Default is "pbcor".
#' @param nboot Number of bootstrap samples. Default 2000.
#' @param alpha Alpha level. Default 0.05. For corfun = pearson, alpha is restricted
#' to 0.05 because the confidence interval adjustments have not been calculated for
#' other alphas.
#' @param alternative Type of test, either "two.sided" (default), "greater" for
#' positive correlations, or "less" for negative correlations.
#' @param null.value Hypothesis to test. Default 0.
#' @param ... Optional parameter to pass to correlation function.
#'
#' @return \itemize{
#' \item \code{estimate.x1y} the correlation coefficient between x1 and y.
#' \item \code{estimate.x2y} the correlation coefficient between x2 and y.
#' \item \code{difference} the difference between correlation coefficients.
#' \item \code{conf.int} the bootstrap confidence interval.
#' \item \code{p.value} the bootstrap p value.
#' \item \code{bootsamples} the bootstrap samples.
#' }
#'
#' @section Note: Modified from function \code{twoDcorR}
#' from Rallfun-v37.txt - see \url{https://github.com/nicebread/WRS/} and
#' \url{http://dornsife.usc.edu/labs/rwilcox/software/}.
#'
#' @references
#' @export
twocorci.ov<-function(x1, x2, y,
                      method="pbcor", nboot=2000, alpha=.05,
                      alternative = "two.sided", null.value=0, ...){
    if(!is.vector(x1) || !is.vector(x2) || !is.vector(y)){
      stop("twocorci.ov: x1, x2, y  must be vectors.")
    }
    if(length(x1)!=length(x2) || length(x1)!=length(y)){
      stop("twocorci.ov: x1, x2, y must have equal lengths.")
    }
    m <- cbind(x1,x2,y)
    m <- m[complete.cases(m), ]
    n <- nrow(m)
    x1 <- m[,1]
    x2 <- m[,2]
    y <- m[,3]
    # correlations
    eval(parse(text=paste("corfun=",method)))
    r.x1y <- corfun(x1,y,...)
    r.x2y <- corfun(x2,y,...)
    # bootstrap
    data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
    data <- listm(t(data))
    bvec <- lapply(data, twocorci.ov.sub, x1, x2, y, corfun, ...) # A 1 by nboot matrix.
    bvec <- matl(bvec)
    bvec <- bvec[1,]-bvec[2,]
    bsort <- sort(bvec)
    ilow <- floor((alpha/2)*nboot+.5)
    ihi <- floor((1-alpha/2)*nboot+.5)
    corci <- c(bsort[ilow], bsort[ihi])
    # p value
    if(alternative == "two.sided"){
      phat <- (sum(bvec < null.value)+.5*sum(bvec==0))/nboot
      sig <- 2 * min(phat, 1 - phat)
    }
    if(alternative == "greater"){
      sig <- 1 - sum(bvec >= null.value)/nboot
    }
    if(alternative == "less"){
      sig <- 1 - sum(bvec <= null.value)/nboot
    }
  list(estimate.x1y=r.x1y, estimate.x2y=r.x2y, difference=r.x1y-r.x2y,
       conf.int=corci, p.value=sig, bootsamples=bvec)
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#' Compare two dependent correlations: non-overlapping case
#'
#' Compute a 1-alpha percentile bootstrap confidence interval for the difference
#' between two correlation coefficients corresponding to two dependent groups,
#' in the non-overlapping case. Compare correlation between x1 and y1 to the correlation between x2 and y2.
#' The default correlation is the percentage bend correlation.
#' This function is inappropriate to make inferences about Pearson's correlations.
#' Missing values are automatically removed.
#'
#' @param x1,y1,x2,y2 Four dependent vectors of the same length.
#' @param method A function that returns a correlation coefficient.
#' Options include "spearman", "pbcor", "wincor". Default is "pbcor".
#' @param nboot Number of bootstrap samples. Default 2000.
#' @param alpha Alpha level. Default 0.05. For corfun = pearson, alpha is restricted
#' to 0.05 because the confidence interval adjustments have not been calculated for
#' other alphas.
#' @param alternative Type of test, either "two.sided" (default), "greater" for
#' positive correlations, or "less" for negative correlations.
#' @param null.value Hypothesis to test. Default 0.
#' @param ... Optional parameter to pass to correlation function.
#'
#' @return \itemize{
#' \item \code{estimate.x1y1} the correlation coefficient between x1 and y1.
#' \item \code{estimate.x2y2} the correlation coefficient between x2 and y2.
#' \item \code{difference} the difference between correlation coefficients.
#' \item \code{conf.int} the bootstrap confidence interval.
#' \item \code{p.value} the bootstrap p value.
#' \item \code{bootsamples} the bootstrap samples.
#' }
#'
#' @section Note: Modified from function \code{twoDNOV}
#' from Rallfun-v37.txt - see \url{https://github.com/nicebread/WRS/} and
#' \url{http://dornsife.usc.edu/labs/rwilcox/software/}.
#'
#' @references
#' @export
twocorci.nov <-function(x1, y1, x2, y2,
                        method="pbcor", nboot=2000, alpha=.05,
                        alternative = "two.sided", null.value=0, ...){
  if(!is.vector(x1) || !is.vector(x2) || !is.vector(y1) || !is.vector(y2)){
    stop("twocorci.nov: x1, x2, y1, y2  must be vectors.")
  }
  if(length(x1)!=length(x2) || length(x1)!=length(y1) || length(y1)!=length(y2)){
    stop("twocorci.nov: x1, x2, y1, y2 must have equal lengths.")
  }
  m <- cbind(x1,y1,x2,y2)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x1 <- m[,1]
  y1 <- m[,2]
  x2 <- m[,3]
  y2 <- m[,4]
  # correlations
  eval(parse(text=paste("corfun=",method)))
  r.x1y1 <- corfun(x1,y1,...)
  r.x2y2 <- corfun(x2,y2,...)
  # bootstrap
  data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
  data <- listm(t(data))
  bvec1 <- matl(lapply(data, corbsub, x1, y1, corfun, ...))
  bvec2 <- matl(lapply(data, corbsub, x2, y2, corfun, ...))
  bvec <- bvec1 - bvec2
  bsort <- sort(bvec)
  ilow <- floor((alpha/2)*nboot+.5)
  ihi <- floor((1-alpha/2)*nboot+.5)
  corci <- c(bsort[ilow], bsort[ihi])
  # p value
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==0))/nboot
    sig <- 2 * min(phat, 1 - phat)
  }
  if(alternative == "greater"){
    sig <- 1 - sum(bvec >= null.value)/nboot
  }
  if(alternative == "less"){
    sig <- 1 - sum(bvec <= null.value)/nboot
  }
  list(estimate.x1y1=r.x1y1, estimate.x2y2=r.x2y2, difference=r.x1y1-r.x2y2,
       conf.int=corci, p.value=sig, bootsamples=bvec)
}
