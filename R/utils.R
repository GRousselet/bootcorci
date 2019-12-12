#' Compute correlation for bootstrap indices
#'
#' This function is used by other functions when computing
#  bootstrap estimates.
#'
#' @param isub A vector of length n. A bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n.
#' @param x,y Two vectors of observations of length n.
#' @param corfun A correlation function which returns a correlation value.
#' @return A bootstrap correlation value.
corbsub <- function(isub, x, y, corfun, ...){
  res <- corfun(x[isub], y[isub], ...)
  res
}

#' Compute correlation for bootstrap indices
#'
#' This function is used by \code{twocorci.ov} when computing
#  bootstrap estimates.
#'
#' @param isub A vector of length n. A bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n.
#' @param x1,x2,y Three vectors of dependent observations of length n.
#' @param corfun A correlation function which returns a correlation value.
#' @return A bootstrap correlation value.
twocorci.ov.sub <- function(isub, x1, x2, y, corfun=pbcor, ...){
  res <- c(corfun(x1[isub], y[isub], ...), corfun(x2[isub], y[isub], ...))
  res
}

#' Adjust quantiles of the correlation bootstrap confidence interval based
#' on sample size. Only necessary for Pearson correlation.
#'
#' Adjustments are available for nboot = 599, so scaling is required for
#' different nboot values. Adjustments work only for alpha = 0.05.
#'
#' @param n Sample size. For the comparison of two independent correlations,
#' n is the sum of the two sample sizes.
#' @param nboot Number of bootstrap samples.
#' @return Indices used to get the 2.5th and 97.5th quantiles of
#' the bootstrap distribution.
adj.corbootci <- function(n, nboot){
  ilow<-15
  ihi<-584
  if(n < 250){
    ilow<-14
    ihi<-585
  }
  if(n < 180){
    ilow<-11
    ihi<-588
  }
  if(n < 80){
    ilow<-8
    ihi<-592
  }
  if(n < 40){
    ilow<-7
    ihi<-593
  }
  round(c(ilow, ihi)*nboot/599)
}

#' Make bivariate normal correlated data
#' @export
mkcord <- function(rho = 0, mu = c(0,0), n = 100){
  cmat <- matrix(c(1, rho, rho, 1), 2, 2)
  out <- MASS::mvrnorm(n = n, mu = mu, Sigma = cmat)
  out
}

#' Store the data in a matrix or data frame in
#' a new R variable having list mode.
#'
#' Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
listm <- function(x){
  if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
  y<-list()
  for(j in 1:ncol(x))y[[j]]<-x[,j]
  y
}

#' Take data in list mode and store it in a matrix
matl <- function(x){
  J=length(x)
  nval=NA
  for(j in 1:J)nval[j]=length(x[[j]])
  temp<-matrix(NA,ncol=J,nrow=max(nval))
  for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
  temp
}
