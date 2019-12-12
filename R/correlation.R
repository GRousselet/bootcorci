#' Compute Spearman's rho
#'
#' Missing values are automatically removed.
#'
#' @param x,y Two vectors of the same length.
#' @param alternative Type of test, either "two.sided", "greater" for positive
#'   correlations, or "less" for negative correlations.
#' @return \itemize{ \item \code{estimate} the correlation coefficient. \item
#' \code{p.value} the p-value of the test. \item \code{statistic} the t
#' statistic of the test. }
#' @export
spearman.test <- function(x, y, alternative="two.sided"){
  if(!is.vector(x) || !is.vector(y)){
    stop("spearman: x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("spearman: the vectors do not have equal lengths.")
  }
  m <- cbind(x,y)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x <- rank(m[,1])
  y <- rank(m[,2])
  # corv <- sum(((x-mean(x))/sqrt(var(x)))*((y-mean(y))/sqrt(var(y)))) / (n-1)
  corv <- cor(x, y)
  test <- corv * sqrt((n - 2)/(1. - corv^2))
  if(alternative == "two.sided"){
    sig <- 2 * (1 - pt(abs(test), n - 2))
  }
  if(alternative == "greater"){
    sig <- 1 - pt(test, n - 2)
  }
  if(alternative == "less"){
    sig <- pt(test, n - 2)
  }
  list(estimate = corv, p.value = sig, statistic = test)
}

#' Compute Pearson's rho
#'
#' Missing values are automatically removed.
#'
#' @param x,y Two vectors of the same length.
#' @return
#' \itemize{
#'   \item \code{estimate} the correlation coefficient.
#'   \item \code{p.value} the p-value of the test.
#'   \item \code{statistic} the t statistic of the test.
#'   }
#' @export
pearson.test <- function(x, y, alternative="two.sided"){
  if(!is.vector(x) || !is.vector(y)){
    stop("pearson: x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("pearson: the vectors do not have equal lengths.")
  }
  m <- cbind(x,y)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x <- m[,1]
  y <- m[,2]
  # corv <- sum(((x-mean(x))/sqrt(var(x)))*((y-mean(y))/sqrt(var(y)))) / (n-1)
  corv <- cor(x, y)
  test <- corv * sqrt((n - 2)/(1 - corv^2))
  if(alternative == "two.sided"){
    sig <- 2 * (1 - pt(abs(test), n - 2))
  }
  if(alternative == "greater"){
      sig <- 1 - pt(test, n - 2)
  }
  if(alternative == "less"){
      sig <- pt(test, n - 2)
  }
  list(estimate = corv, p.value = sig, statistic = test)
}

#' Compute the percentage bend correlation between x and y.
#'
#' @param x,y Two vectors of the same length.
#' @param beta The bending constant (default 0.2).
#' @return
#' \itemize{
#'   \item \code{estimate} the correlation coefficient.
#'   \item \code{p.value} the p-value of the test.
#'   \item \code{statistic} the t statistic of the test.
#'   }
#' @export
pbcor.test <- function(x, y, beta=.2, alternative="two.sided"){
  if(!is.vector(x) || !is.vector(y)){
    stop("pbcor: x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("pbcor: the vectors do not have equal lengths.")
  }
  m <- cbind(x,y)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x <- m[,1]
  y <- m[,2]
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
  test <- corv*sqrt((n - 2)/(1 - corv^2))
  if(alternative == "two.sided"){
    sig <- 2 * (1 - pt(abs(test), n - 2))
  }
  if(alternative == "greater"){
    sig <- 1 - pt(test, n - 2)
  }
  if(alternative == "less"){
    sig <- pt(test, n - 2)
  }
  list(estimate=corv, statistic=test, p.value=sig)
}

#' Compute the Winsorized correlation between x and y.
#'
#' @param x,y Two vectors of the same length.
#' @param tr The amount of Winsorization (default 0.2).
#' @return
#' \itemize{
#'   \item \code{estimate} the correlation coefficient.
#'   \item \code{p.value} the p-value of the test.
#'   \item \code{statistic} the t statistic of the test.
#'   }
#' @export
wincor.test <- function(x, y, tr=0.2, alternative="two.sided"){
  if(!is.vector(x) || !is.vector(y)){
    stop("wincor: x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("wincor: the vectors do not have equal lengths.")
  }
  m <- cbind(x,y)
  m <- m[complete.cases(m), ]
  n <- nrow(m)
  x <- m[,1]
  y <- m[,2]
  g <- floor(tr*n)
  xvec <- winval(x,tr)
  yvec <- winval(y,tr)
  corv <- cor(xvec,yvec)
  test <- corv*sqrt((n-2)/(1.-corv^2))
  if(alternative == "two.sided"){
    sig <- 2 * (1 - pt(abs(test), n-2*g-2))
  }
  if(alternative == "greater"){
    sig <- 1 - pt(test, n-2*g-2)
  }
  if(alternative == "less"){
    sig <- pt(test, n-2*g-2)
  }
  list(estimate=corv, statistic=test, p.value=sig)
}

#' Winsorize the data in the vector x.
#'
#' Called by \code{wincor}.
#'
#' @param x A vector.
#' @param tr The amount of Winsorization, between 0 and 1 (default 0.2).
#'
#' @export
winval <- function(x, tr=0.2){
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr*n)+1
  itop <- length(x)-ibot+1
  xbot <- y[ibot]
  xtop <- y[itop]
  winval <- ifelse(x<=xbot,xbot,x)
  winval <- ifelse(winval>=xtop,xtop,winval)
  winval
}

#' Compute the one-step percentage bend measure of location
pbos <- function(x, beta=.2){
  n <- length(x)
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  psi <- (x-median(x))/omhatx
  i1 <- length(psi[psi<(-1)])
  i2 <- length(psi[psi>1])
  sx <- ifelse(psi<(-1),0,x)
  sx <- ifelse(psi>1,0,sx)
  pbos <- (sum(sx) + omhatx*(i2-i1)) / (n-i1-i2)
  pbos
}
