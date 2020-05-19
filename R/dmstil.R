#' This function return the quasi density function for mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param lambda skewing matrix with k rows.
#' @param delta location vector of size k.
#' @param Ainv lower triangular matrix, where t(Ainv) * Ainv is the precision matrix.
#' @param nu degree of freedom (>0).
#' @param u (Optional) a m x k matrix, a set of samples generated from standard k-dimensional multivariate t distribution with degree of freedom nu.
#' @param log.p a logical value. If TRUE, return the probability density function in logarithmic scale. By default FALSE.
#' @param control a list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6.}
##' }
#' @return return a numeric vector of length n.
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # k <- ncol(RiverFlow)
#' # lambda <- diag(k)
#' # delta <- rep(0,k)
#' # Ainv <- diag(k)
#' # nu <- 2
#' # dmstil(as.matrix(log(RiverFlow)), lambda, delta, Ainv, nu)
dmstil <- function(x, lambda, delta, Ainv, nu, u, log.p = FALSE, control = list()) {
  .check.control(control)
  if (!"numLikSample" %in% names(control)) control$numLikSample <- 1e6
  numLikSample <- control$numLikSample
  
  k <- ncol(x)
  .check.mstil.param(k, lambda, delta, Ainv, nu)
  
  if (missing(u)) u <- mvtnorm::rmvt(numLikSample, delta = rep(0, k), sigma = diag(k), df = nu)
  
  z <- t((t(x) - delta)) %*% t(Ainv)
  
  Gz <- stats::plogis(z %*% lambda, log.p = TRUE)
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  
  logDensity <- rowSums(Gz) + .dmvt2(x, delta = delta, Ainv = Ainv, nu = nu, log.p = TRUE) - log(mean(exp(rowSums(Gu))))
  
  if (log.p) {
    return(logDensity)
  } else {
    return(exp(logDensity))
  }
}
