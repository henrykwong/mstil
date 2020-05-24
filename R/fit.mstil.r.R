#' This function finds the maximum likelihood estimations for restricted mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of inital parameters, contains lambda, delta, Ainv, and nu.
#' @param control list of control variables, see 'details'.
#' #@details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{lambdaPenalty}{a positive value, represents the L2 penalty coefficient for lambda. By default 0.}
##'  \item{ainvPenalty}{a positive value, represents the L2 penalty coefficient for Ainv. By default 0.}
##'  \item{maxitOptimR}{a positive integer, represents the maximum number of iterations allowed in optim. By default 1e3.}
##' }
#' @return a list with components:
#' \item{lambda}{the value of fitted lambda.}
#' \item{delta}{the value of fitted delta.}
#' \item{Ainv}{the value of fitted Ainv.}
#' \item{nu}{the value of fitted nu.}
#' \item{logLik}{the value of fitted log-Likelihood.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.mstil.r(as.matrix(log(RiverFlow)))
fit.mstil.r <- function(x, param = NULL, control = list()) {
  
  .check.control(control)
  
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 0
  if (!"maxitOptimR" %in% names(control)) control$maxitOptimR <- 1e3
  if (!"ainvPenalty" %in% names(control)) control$ainvPenalty <- 0
  ainvPenalty <- nrow(x) * control$ainvPenalty
  lambdaPenalty <- nrow(x) * control$lambdaPenalty
  maxitOptimR <- control$maxitOptimR
  k <- ncol(x)
  n <- nrow(x)
  if (missing(param) | is.null(param)) param <- .default.init.param.method.t(x)
  .check.mstil.r.param(k, param$lambda, param$delta, param$Ainv, param$nu)
  
  lambda <- param$lambda
  delta <- param$delta
  Ainv <- param$Ainv
  nu <- param$nu
  
  lik <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    lnu <- param[length(param)]
    nu <- exp(lnu)
    res <- sum(dmstil.r(x, lambda, delta, Ainv, nu, log.p = TRUE)) - lambdaPenalty * sum(lambda ^ 2) - ainvPenalty * sum(Ainv ^ 2)
    return(res)
  }
  grad <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    lnu <- param[length(param)]
    nu <- exp(lnu)
    grad <- .mstil.r.grad(x, lambda, delta, Ainv, nu)
    grad$dlambda <- grad$dlambda - 2 * lambdaPenalty * lambda
    grad$dAinv <- grad$dAinv - 2 * ainvPenalty * Ainv
    gr <- c(diag(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)], grad$dlnu)
    
    return(gr)
  }
  param0 <- c(diag(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)], log(nu))
  
  res <- stats::optim(param0, lik, grad, method = "BFGS", control = c(fnscale = -1, maxit = maxitOptimR))
  param1 <- res$par
  lambda1 <- diag(param1[1:k])
  delta1 <- param1[1:k + (k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (k) + 1):(length(param1) - 1)]
  lnu1 <- param1[length(param1)]
  nu1 <- exp(lnu1)
  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu1, logLik = res$value))
}
