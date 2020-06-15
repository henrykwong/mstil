#' This function estimate the value of log likelihood function
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param lambda skewing matrix with k rows.
#' @param delta location vector of size k.
#' @param Ainv lower triangular matrix, where t(Ainv) * Ainv is the precision matrix.
#' @param nu degree of freedom (>0).
#' @param u (Optional) a m x k matrix, a set of samples generated from standard k-dimensional multivariate t distribution with degree of freedom nu.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{a positive integer, represents the number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{a value between 0.5 and 1, represents the 1 sided confidence level of the log-likelihood to be calculated. By default 0.95.}
##' }
#' @return a list with components:
#' \item{logLikLower}{the lower bound of the estimated log-likelihood function.}
#' \item{logLik}{the estimate of the log-likelihood function.}
#' \item{logLikUpper}{the upper bound of the estimated log-likelihood function.}
#' @export
mstil.logLik <- function(x, lambda, delta, Ainv, nu, u, control = list()) {
  n <- nrow(x)
  k <- ncol(x)
  .check.control(control)
  .check.mstil.param(k, lambda, delta, Ainv, nu)
  
  if (!"numLikSample" %in% names(control)) control$numLikSample <- 1e6
  if (!"conLevel" %in% names(control)) control$conLevel <- 0.95
  numLikSample <- control$numLikSample
  conLevel <- control$conLevel
  

  
  if (missing(u)) u <- mvtnorm::rmvt(numLikSample, delta = rep(0, k), sigma = diag(k), df = nu)
  
  z <- t((t(x) - delta)) %*% t(Ainv)
  Gz <- stats::plogis(z %*% lambda, log.p = TRUE)
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  
  expGu <- exp(rowSums(Gu))
  EGu <- mean(expGu)
  spread <- abs(stats::qnorm((1 - conLevel))) * stats::sd(expGu) / sqrt(nrow(u))
  
  res <- sum(rowSums(Gz) + .dmvt2(x, delta = delta, Ainv = Ainv, nu = nu, log.p = TRUE))
  
  return(list(logLikLowerBound = res - n * log(EGu + spread), logLik = res - n * log(EGu)))
}