#' This function finds the maximum likelihood estimations for restricted mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of inital parameters, contains lambda, delta, Ainv, and nu.
#' @param control list of control variables, see 'details'.
#' #@details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{lambdaPenalty}{a positive value, represents the L1 penalty coefficient for lambda. By default 0.}
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
fit.mstil.r <- function(x, param, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 0
  if (!"maxitOptimR" %in% names(control)) control$maxitOptimR <- 1e3
  lambdaPenalty <- control$lambdaPenalty
  maxitOptimR <- control$maxitOptimR
  k <- ncol(x)
  n <- nrow(x)
  
  if (missing(param)) {
    param <- list(
      lambda = 0 * diag(stats::runif(k, -0.1, 0.1)),
      delta = colMeans(x),
      Ainv = tryCatch(t(solve(chol(stats::cov(x)))),
                      error = function(e) 1 / sqrt(stats::cov(x)),
                      warning = function(w) 1 / sqrt(stats::cov(x))
      ),
      nu = 10
    )
    param$Ainv[upper.tri(param$Ainv, diag = FALSE)] <- 0
  }
  
  lambda <- param$lambda
  delta <- param$delta
  Ainv <- param$Ainv
  nu <- param$nu
  
  if (any((lambda * diag(k)) != lambda)) stop("lambda is not a diagonal matrix!")
  if (nrow(lambda) != k) stop("lambda is non-conformable!")
  if (length(delta) != k) stop("length of delta is not equal to dimension of x!")
  if (any(dim(Ainv) != k) | any(Ainv[upper.tri(Ainv)] != 0)) stop("Ainv is non-conformable or is not an lower triangular matrix!")
  if (nu <= 0) stop("nu is non-positive!")
  
  lik <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    lnu <- param[length(param)]
    nu <- exp(lnu)
    res <- sum(dmstil.r(x, lambda, delta, Ainv, nu, log.p = TRUE)) - lambdaPenalty * sum(abs(lambda))
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
    grad$dlambda <- grad$dlambda - lambdaPenalty * sign(lambda)
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
