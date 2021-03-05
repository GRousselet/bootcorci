# Simplified versions of the correlation functions to use with the bootstrap functions.
# Only the correlation coefficients are returned and no checks are done.

#' Compute Kendall's tau
#'
#' Simplified version of \code{kendall.test} which returns only the correlation coefficient.
#' Called by the bootstrap functions for faster execution.
#'
#' @param x,y Two vectors of the same length.
#' @return The correlation coefficient.
#' @export
kendall <- function(x,y){
  n <- length(x)
  xdif <- outer(x,x,FUN="-")
  ydif <- outer(y,y,FUN="-")
  tv <- sign(xdif)*sign(ydif)
  corv <- sum(tv)/(n*(n-1))
  corv
}

#' Compute Spearman's rho
#'
#' Simplified version of \code{spearman.test} which returns only the correlation coefficient.
#' Called by the bootstrap functions for faster execution.
#'
#' @param x,y Two vectors of the same length.
#' @return The correlation coefficient.
#' @export
spearman <- function(x,y){
  n <- length(x)
  x <- rank(x)
  y <- rank(y)
  # corv <- sum(((x-mean(x))/sqrt(var(x)))*((y-mean(y))/sqrt(var(y)))) / (n-1)
  corv <- cor(x, y)
  corv
}

#' Compute Pearson's rho
#'
#' Simplified version of \code{pearson.test} which returns only the correlation coefficient.
#' Called by the bootstrap functions for faster execution.
#'
#' @param x,y Two vectors of the same length.
#' @return The correlation coefficient.
#' @export
pearson <- function(x,y){
  n <- length(x)
  # corv <- sum(((x-mean(x))/sqrt(var(x)))*((y-mean(y))/sqrt(var(y)))) / (n-1)
  corv <- cor(x, y)
  corv
}

#' Compute the percentage bend correlation between x and y.
#'
#' Simplified version of \code{pbcor.test} which returns only the correlation coefficient.
#' Called by the bootstrap functions for faster execution.
#'
#' @param x,y Two vectors of the same length.
#' @param beta The bending constant (default 0.2).
#' @return The correlation coefficient.
#' @export
pbcor <- function(x, y, beta=.2){
  n <- length(x)
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  temp <- sort(abs(y-median(y)))
  omhaty <- temp[floor((1-beta)*n)]
  a <- (x-pbos(x,beta))/omhatx
  b <- (y-pbos(y,beta))/omhaty
  a <- ifelse(a<=-1,-1,a)
  a <- ifelse(a>=1,1,a)
  b <- ifelse(b<=-1,-1,b)
  b <- ifelse(b>=1,1,b)
  corv <- sum(a*b)/sqrt(sum(a^2)*sum(b^2))
  corv
}

#' Compute the Winsorized correlation between x and y.
#'
#' Simplified version of \code{wincor.test} which returns only the correlation coefficient.
#' Called by the bootstrap functions for faster execution.
#'
#' @param x,y Two vectors of the same length.
#' @param tr The amount of Winsorization (default 0.2).
#' @return The correlation coefficient.
#'
#' @export
wincor <- function(x, y, tr=0.2){
  n <- length(x)
  g <- floor(tr*n)
  x <- winval(x,tr)
  y <- winval(y,tr)
  corv <- cor(x, y)
  corv
}

