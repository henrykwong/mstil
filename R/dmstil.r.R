#' This function returns to likelihood function of restricted mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param lambda skewing matrix with k rows.
#' @param delta location vector of size k.
#' @param Ainv lower triangular matrix, where t(Ainv) * Ainv is the precision matrix.
#' @param nu degree of freedom (>0).
#' @param log.p a logical value. If TRUE, return the probability density function in logarithmic scale. By default FALSE.
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
#' # dmstil.r(as.matrix(log(RiverFlow)), lambda, delta, Ainv, nu)
dmstil.r <- function(x, lambda, delta, Ainv, nu, log.p = FALSE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (any(diag(diag(lambda)) != lambda)) {
    warning("lambda must be a diagonal matrix!")
    return(rep(NaN, nrow(x)))
  }
  .check.mstil.r.param(ncol(x), lambda, delta, Ainv, nu)
  k <- ncol(x)
  z <- t((t(x) - delta)) %*% t(Ainv)
  Gz <- stats::plogis(z %*% lambda, log.p = TRUE)
  res <- rowSums(Gz) + .dmvt2(x, delta = delta, Ainv = Ainv, nu = nu, log.p = TRUE) + k * log(2)
  if (log.p) {
    return(res)
  } else {
    return(exp(res))
  }
}
