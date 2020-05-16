#' This function returns the log likelihood, ICL, AIC, BIC, and group memberships for fmmstil
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of lists of parameters, contains list of omega, list of lambda, list delta, list of Ainv, and list of nu.
#' @param u (Optional) list of K m x k matrices, each matrix contain samples generated from standard k-dimensional multivariate t distribution with degree of freedom of the K-th cluster.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6.}
##' }
#' @return a list with components:
#' \item{logLik}{the value of log-likelihood.}
#' \item{ICL}{the value of integrated completed log-likelihood.}
#' \item{AIC}{the value of Akaike information criterion.}
#' \item{BIC}{the value of Bayesian information criterion.}
#' \item{clust}{the estimated group membership.}
#' @export
fmmstil.fitness <- function(x, param, u, control = list()) {
  weight <- .fmmstil.weight(x, param$omega, param$lambda,param$delta, param$Ainv, param$nu, u = u, control = control)
  K <- length(param$omega)
  k <- ncol(x)
  m <- k * ( k + 1) + 1 + k
  n <- nrow(x)
  logLik <- sum(log(rowSums(weight)))
  weight <- weight / rowSums(weight)
  guess <- apply(weight, 1, which.max)
  nK <- table(guess)
  mK <- m * K + K - 1
  ICL <- logLik +
    mK / 2 * log(n) +
    lgamma(K / 2) + sum(log(apply(weight, 1, max)))
  sum(lgamma(nK + 0.5)) -
    K * lgamma(0.5) -
    lgamma(n + K / 2)
  if (any(nK < 10) || length(nK) < K) ICL <- -Inf
  return( list(ICL = ICL, clust = guess, BIC = log(n) * mK - 2 * logLik, AIC <- 2 * mK - 2 * logLik))
}
