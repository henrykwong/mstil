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
  if (!"numLikSample" %in% names(control)) control$numLikSample <- 1e6
  numLikSample <- control$numLikSample
  k <- ncol(x)
  if (missing(u)) u <- mvtnorm::rmvt(numLikSample, delta = rep(0, k), sigma = diag(k), df = nu)
  if (nrow(lambda) != k) stop("number of rows of lambda is not equal to dimension of x!")
  if (length(delta) != k) stop("length of delta is not equal to dimension of x!")
  if (ncol(Ainv) != nrow(Ainv) | length(Ainv) != k^2 | any(Ainv[upper.tri(Ainv)] != 0)) stop("Ainv is of wrong size or is not an lower triangular matrix!")
  if (nu < 0) stop("nu must be positive!")

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
  k <- ncol(x)
  if (nrow(lambda) != k) {
    warning("number of rows of lambda is not equal to dimension of x!")
    return(rep(NaN, nrow(x)))
  }
  if (length(delta) != k) {
    warning("length of delta is not equal to dimension of x!")
    return(rep(NaN, nrow(x)))
  }
  if (ncol(Ainv) != nrow(Ainv) | length(Ainv) != k^2 | any(Ainv[upper.tri(Ainv)] != 0)) {
    warning("Ainv is of wrong size or is not an lower triangular matrix!")
    return(rep(NaN, nrow(x)))
  }
  if (nu < 0) {
    warning("nu must be positive!")
    return(rep(NaN, nrow(x)))
  }

  z <- t((t(x) - delta)) %*% t(Ainv)
  Gz <- stats::plogis(z %*% lambda, log.p = TRUE)
  res <- rowSums(Gz) + .dmvt2(x, delta = delta, Ainv = Ainv, nu = nu, log.p = TRUE) + k * log(2)
  if (log.p) {
    return(res)
  } else {
    return(exp(res))
  }
}



#' This function estimate the value of log likelihood function
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param lambda skewing matrix with k rows.
#' @param delta location vector of size k.
#' @param Ainv lower triangular matrix, where t(Ainv) * Ainv is the precision matrix.
#' @param nu degree of freedom (>0).
#' @param u (Optional) a m x k matrix, a set of samples generated from standard k-dimensional multivariate t distribution with degree of freedom nu.
#' @param control a list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{numeric between 0.5 and 1. Confidence level of the log-likelihood. By default 0.95.}
##' }
#' @return a list with components:
#' \item{logLikLower}{the lower bound of the estimated log-likelihood function.}
#' \item{logLik}{the estimate of the log-likelihood function.}
#' \item{logLikUpper}{the upper bound of the estimated log-likelihood function.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # k <- ncol(RiverFlow)
#' # lambda <- diag(k)
#' # delta <- rep(0,k)
#' # Ainv <- diag(k)
#' # nu <- 2
#' # mstil.logLik(as.matrix(log(RiverFlow)), lambda, delta, Ainv, nu)
mstil.logLik <- function(x, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numLikSample" %in% names(control)) control$numLikSample <- 1e6
  if (!"conLevel" %in% names(control)) control$conLevel <- 0.95
  numLikSample <- control$numLikSample
  conLevel <- control$conLevel
  k <- ncol(x)
  if (nrow(lambda) != k) stop("number of rows of lambda is not equal to dimension of x!")
  if (length(delta) != k) stop("length of delta is not equal to dimension of x!")
  if (ncol(Ainv) != nrow(Ainv) | length(Ainv) != k^2 | any(Ainv[upper.tri(Ainv)] != 0)) stop("Ainv is of wrong size or is not an lower triangular matrix!")
  if (nu < 0) stop("nu must be positive!")
  if (conLevel >= 1 | conLevel <= 0.5) stop("conLevel must be between 0.5 and 1!")
  if (is.data.frame(x)) x <- as.matrix(x)
  n <- nrow(x)
  z <- t((t(x) - delta)) %*% t(Ainv)

  if (missing(u)) u <- mvtnorm::rmvt(numLikSample, delta = rep(0, k), sigma = diag(k), df = nu)

  Gz <- stats::plogis(z %*% lambda, log.p = TRUE)
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)

  expGu <- exp(rowSums(Gu))
  EGu <- mean(expGu)
  spread <- abs(stats::qnorm((1 - conLevel) / 2)) * stats::sd(expGu) / sqrt(nrow(u))

  res <- sum(rowSums(Gz) + .dmvt2(x, delta = delta, Ainv = Ainv, nu = nu, log.p = TRUE))
  return(list(logLikLower = res - n * log(EGu + spread), logLik = res - n * log(EGu), logLikUpper = res - n * log(EGu - spread)))
}



#' This function finds the maximum penalised quasi likelihood estiamtes for mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of inital parameters, containing lambda, delta, Ainv, and nu.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{numeric between 0.5 and 1. Confidence level of the log-likelihood. By default 0.95.}
##'  \item{cvgCon}{non-negative interger. The algorithm stops when the estimated maximum penalised log-likelihood is not improved in cvgCon iterations. By default 5.}
##'  \item{lambdaPenalty}{L1 penalty coefficient for lambda. By default 0.}
##'  \item{maxit}{maximum number iterations. By default 1000.}
##'  \item{maxitOptim}{maximum number of iterations in optim. By default 10.}
##'  \item{numGradSample}{number of samples used to estimate the gradient. By default 1e4.}
##'  \item{finDiffStep}{step size to estimate gradient w.r.t. nu. By default 1e-5.}
##'  \item{stepSgd}{step size to be used in the stochastic gradient step. By default 1e-2.}
##'  \item{iterSgd}{number of iterations to be used in the stochastic gradient step. By default 1e2.}
##'  \item{dimRateSgd}{dimishing rate for step size to be used in the stochastic gradient step. By default 1e-2.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector of the estimated log-likelihood function after each itereation.}
#' \item{par}{a list of list of fitted parameters after each iteration. Within each list, "lambda" is the fitted skewing matrix, "delta" is the location parameter, "Ainv" is the reparameterised scale parameter, "nu" is the degree of freedom.}
#' \item{logLikLower}{a numeric vector of the lower bound of the estimated log-likelihood function after each iteration.}
#' \item{logLikUpper}{a numeric vector of the upper bound of the estimated log-likelihood function after each iteration.}
#' \item{time}{a non-negative numeric vector, records the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.mstil(as.matrix(log(RiverFlow)))
fit.mstil <- function(x, param, show.progress = TRUE, control = list()) {
  if (!"maxit" %in% names(control)) control$maxit <- 1e3
  if (!"cvgCon" %in% names(control)) control$cvgCon <- 5

  maxit <- control$maxit
  cvgCon <- control$cvgCon

  k <- ncol(x)
  if (missing(param)) param <- fit.mstil.r(x, control = control)
  if (nrow(param$lambda) != k) stop("Number of rows of lambda is not equal to dimension of x!")
  if (length(param$delta) != k) stop("Length of delta is not equal to dimension of x!")
  if (ncol(param$Ainv) != nrow(param$Ainv) | length(param$Ainv) != k^2 | any(param$Ainv[upper.tri(param$Ainv)] != 0)) stop("Ainv is of wrong size or is not an lower triangular matrix!")
  if (param$nu < 0) stop("nu must be positive!")


  res <- list(lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  
  startTime <- Sys.time()
  
  lik <- mstil.logLik(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
  likRec <- lik$logLik
  likLowerRec <- lik$logLikLower
  likUpperRec <- lik$logLikUpper

  timeRec <- 0

  resRec <- list()
  resRec[[1]] <- res


  for (i in 2:(maxit + 1)) {
    res <- .fit.mstil.1(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    res$nu <- .fit.mstil.2(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    lik <- mstil.logLik(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    resRec[[i]] <- res
    likRec <- c(likRec, lik$logLik)
    likLowerRec <- c(likLowerRec, lik$logLikLower)
    likUpperRec <- c(likUpperRec, lik$logLikUpper)
    timeRec <- c(timeRec, as.numeric(Sys.time() - startTime, units = "secs"))
    if (show.progress) {
      cat("\r", "Iteration : ", (i - 1), "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }

    if (i > cvgCon) {
      if (all(likRec[i - cvgCon] > likRec[(i - cvgCon + 1):i])) {
        if (show.progress) cat("\n", "Converged!")
        return(list(logLik = likRec, par = resRec, logLikLower = likLowerRec, logLikUpper = likUpperRec, time = timeRec))
        break
      }
    }
  }
  if (show.progress) cat("\n", "Maximum number of iterations reached!")
  return(list(logLik = likRec, par = resRec, logLik_lower = likLowerRec, logLik_upper = likUpperRec, time = timeRec))
}



#' This function finds the maximum likelihood estimations for restricted mstil.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param param list of inital parameters, containing lambda, delta, Ainv, and nu.
#' @param control list of control variables, see 'details'.
#' #@details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{lambdaPenalty}{L1 penalty coefficient for lambda. By default 0.}
##'  \item{maxitOptimR}{maximum number of iterations in optim.. By default 1e3.}
##' }
#' @return a list with components:
#' \item{lambda}{fitted lambda.}
#' \item{delta}{fitted delta.}
#' \item{Ainv}{fitted Ainv.}
#' \item{nu}{fitted nu.}
#' \item{logLik}{value of fitted log-Likelihood.}
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

  if (any((lambda * diag(k)) != lambda)) stop("lambda must be a diagonal matrix!")
  if (nrow(lambda) != k) stop("number of rows of lambda is not equal to dimension of x!")
  if (length(delta) != k) stop("length of delta is not equal to dimension of x!")
  if (ncol(Ainv) != nrow(Ainv) | length(Ainv) != k^2 | any(Ainv[upper.tri(Ainv)] != 0)) stop("Ainv is of wrong size or is not an lower triangular matrix!")
  if (nu < 0) stop("nu must be positive!")

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



#' This function finds the penalised quasi likelihood estiamtes for finite mixture of mstil via EM.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param K number of clusters.
#' @param param list of lists of inital parameters, containing omega, lambda, delta, Ainv, and nu.
#' @param init.cluster (optional) initial clusters used to find initial parameters.
#' @param init.param.method (optional) method to obtain initial parameters. It needs to be a function of x and K, and return a list of list of parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{numeric between 0.5 and 1. Confidence level of the log-likelihood. By default 0.95.}
##'  \item{cvgCon}{non-negative interger. The algorithm stops when the estimated maximum penalised log-likelihood is not improved in cvgCon iterations. By default 5.}
##'  \item{lambdaPenalty}{L1 penalty coefficient for lambda. By default 0.}
##'  \item{maxit}{maximum number of EM iterations. By default 1000.}
##'  \item{maxitOptim}{maximum number of iterations in optim within each M-step. By default 10.}
##'  \item{numGradSample}{number of samples used to estimate the gradient. By default 1e4.}
##'  \item{finDiffStep}{step size to estimate gradient w.r.t. nu. By default 1e-5.}
##'  \item{stepSgd}{step size to be used in the stochastic gradient step. By default 1e-2.}
##'  \item{iterSgd}{number of iterations to be used in the stochastic gradient step. By default 1e2.}
##'  \item{dimRateSgd}{dimishing rate for step size to be used in the stochastic gradient step. By default 1e-2.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector recording the estimated log-likelihood after each itereation.}
#' \item{par}{a list of list of list of fitted parameters after each iteration.}
#' \item{time}{a vectorrecording the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.fmmstil(as.matrix(log(RiverFlow)), 2)
fit.fmmstil <- function(x, K, param, init.cluster, init.param.method, show.progress = TRUE, control = list()) {
  if (!"maxit" %in% names(control)) control$maxit <- 1e3
  if (!"cvgCon" %in% names(control)) control$cvgCon <- 5
  maxit <- control$maxit
  cvgCon <- control$cvgCon
  n <- nrow(x)
  k <- ncol(x)
  res <- list()
  if (missing(param)) {
    param <- list(omega = list(), lambda = list(), delta = list(), Ainv = list(), nu = list())
    if (missing(init.param.method)) init.param.method <- .default.init.param.method
    if (missing(init.cluster)) init.cluster <- .default.init.cluster.method(x, K)
    param$omega <- as.list(table(init.cluster) / n)
    for (i in 1:K) {
      initFit <- init.param.method(x[which(init.cluster == unique(init.cluster)[i]), ])
      param$lambda[[i]] <- initFit$lambda
      param$delta[[i]] <- initFit$delta
      param$Ainv[[i]] <- initFit$Ainv
      param$nu[[i]] <- initFit$nu
    }
  }
  
  res <- list(omega = param$omega, lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  
  startTime <- Sys.time()
  likRec <- c()
  resRec <- list()
  timeRec <- c()
  for (i in 1:maxit) {
    w_ <- .fmmstil.weight(x, res$omega, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    d <- rowSums(w_)
    logLik <- sum(log(d))
    w <- w_ / d
    res$omega <- as.list(colSums(w) / n)
    likRec <- c(likRec, logLik)
    resRec[[i]] <- res
    timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))
    if (i > cvgCon) {
      if (all(likRec[i - cvgCon] > likRec[(i - cvgCon + 1):i])) {
        if (show.progress) cat("\n", "converged!")
        return(list(logLik = likRec, par = resRec, time = timeRec))
        break
      }
    }
    if (show.progress) {
      cat("\r", "Iteration : ", i, "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    for (j in 1:K) {
      res1 <- .fit.mstil.1.weighted(x, w[, j], lambda = res$lambda[[j]], delta = res$delta[[j]], Ainv = res$Ainv[[j]], nu = res$nu[[j]], control = control)
      res$lambda[[j]] <- res1$lambda
      res$delta[[j]] <- res1$delta
      res$Ainv[[j]] <- res1$Ainv
      res$nu[[j]] <- .fit.mstil.2.weighted(x, w[, j], lambda = res$lambda[[j]], delta = res$delta[[j]], Ainv = res$Ainv[[j]], nu = res$nu[[j]], control = control)
    }
  }

  w_ <- .fmmstil.weight(x, res$omega, res$lambda, res$delta, res$Ainv, res$nu, control = control)
  d <- rowSums(w_)
  logLik <- sum(log(d))
  w <- w_ / d
  res$omega <- as.list(colSums(w) / n)
  likRec <- c(likRec, logLik)
  resRec[[length(likRec)]] <- res
  timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))

  if (show.progress) {
    cat("\r", "Iteration : ", length(likRec), "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    cat("\n", "Maximum number of iteration reached!")
  }
  return(list(logLik = likRec, par = resRec, time = timeRec))
}



#' This function finds the maximum likelihood estiamtes for finite mixture of restricted mstil via EM.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param K number of clusters.
#' @param param list of lists of inital parameters, containing omega, lambda, delta, Ainv, and nu.
#' @param init.cluster (optional) initial clusters used to find initial parameters.
#' @param init.param.method (optional) method to obtain initial parameters. It needs to be a function of x that returns a list of list of parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{cvgConR}{tolerance level. The algorithm stops when increase of log-likelihood drops below this level. }
##'  \item{lambdaPenalty}{L1 penalty coefficient for lambda. By default 0.}
##'  \item{maxitR}{maximum number of EM iterations. By default 1000.}
##'  \item{maxitOptimR}{maximum number of iterations in optim within each M-step. By default 1e3.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector recording the estimated log-likelihood after each itereation.}
#' \item{par}{a list of list of list of fitted parameters after each iteration.}
#' \item{time}{a vectorrecording the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.fmmstil.r(as.matrix(log(RiverFlow)), 2)
fit.fmmstil.r <- function(x, K, param, init.cluster, init.param.method, show.progress = TRUE, control = list()) {
  if (!"maxitR" %in% names(control)) control$maxitR <- 1e3
  if (!"cvgConR" %in% names(control)) control$cvgConR <- 1e-2
  maxitR <- control$maxitR
  cvgConR <- control$cvgConR
  n <- nrow(x)
  k <- ncol(x)
  res <- list()

  if (missing(param)) {
    param <- list(omega = list(), lambda = list(), delta = list(), Ainv = list(), nu = list())
    if (missing(init.param.method)) init.param.method <- .default.init.param.method
    if (missing(init.cluster)) init.cluster <- .default.init.cluster.method(x, K)
    param$omega <- as.list(table(init.cluster) / n)
    for (i in 1:K) {
      initFit <- init.param.method(x[which(init.cluster == unique(init.cluster)[i]), ])
      param$lambda[[i]] <- initFit$lambda
      param$delta[[i]] <- initFit$delta
      param$Ainv[[i]] <- initFit$Ainv
      param$nu[[i]] <- initFit$nu
    }
  }

  res <- list(omega = param$omega, lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  startTime <- Sys.time()
  likRec <- c()
  resRec <- list()
  timeRec <- c()
  for (i in 1:maxitR) {
    w_ <- .fmmstil.r.weight(x, res$omega, res$lambda, res$delta, res$Ainv, res$nu)
    d <- rowSums(w_)
    logLik <- sum(log(d))
    w <- w_ / d
    res$omega <- as.list(colSums(w) / n)
    likRec <- c(likRec, logLik)
    resRec[[i]] <- res
    timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))
    if ((likRec[i] - max(-Inf, likRec[i - 1])) < cvgConR) {
      if (show.progress) cat("\n", "Converged!")
      return(list(logLik = likRec, par = resRec, time = timeRec))
      break
    }


    if (show.progress) {
      cat("\r", "Iteration : ", i, "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    for (j in 1:K) {
      res1 <- .fit.fmmstil.r.weighted(x, w[, j], lambda = res$lambda[[j]], delta = res$delta[[j]], Ainv = res$Ainv[[j]], nu = res$nu[[j]], control = control)
      res$lambda[[j]] <- res1$lambda
      res$delta[[j]] <- res1$delta
      res$Ainv[[j]] <- res1$Ainv
      res$nu[[j]] <- res1$nu
    }
  }

  w_ <- .fmmstil.r.weight(x, res$omega, res$lambda, res$delta, res$Ainv, res$nu)
  d <- rowSums(w_)
  logLik <- sum(log(d))
  w <- w_ / d
  res$omega <- as.list(colSums(w) / n)
  likRec <- c(likRec, logLik)
  resRec[[length(likRec)]] <- res
  timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))

  if (show.progress) {
    cat("\r", "Iteration : ", length(likRec), "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    cat("\n", "Maximum number of iteration reached!")
  }
  return(list(logLik = likRec, par = resRec, time = timeRec))
}



#' Clustering using mstil given a fixed number of clusters.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param K positive integer, number of cluster.
#' @param numTrial a positive integer, number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{restricted}{a list containing details of the best fitted fmmstil.r.}
#' \item{unrestricted}{a list containing details of the best fitted fmmstil.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.K(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.K <- function(x, K, numTrial, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgConR" %in% names(control)) control$cvgConR <- 1e-1
    
    
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  
  resRec <- list()

  maxICL <- -Inf
  for (i in 1:numTrial) {
    if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Trial : ", i, " of ", numTrial)

    init.cluster <- init.cluster.method(x, K)

    res1 <- tryCatch(fit.fmmstil.r(x, K, init.cluster = init.cluster, init.param.method = init.param.method, show.progress = FALSE, control = control),
      error = function(e) NA,
      warning = function(w) NA
    )
    if (is.list(res1)) {
      par <- res1$par[[which.max(res1$logLik)]]
      fitness1 <- fmmstil.r.fitness(x, par)
      res1$ICL <- fitness1$ICL
      res1$clust <- fitness1$clust
      res1$AIC <- fitness1$AIC
      res1$BIC <- fitness1$BIC
      
      if (res1$ICL > maxICL) {
        res1Best <- res1
        maxICL <- res1$ICL
        guess1Best <- res1$clust
      }
      
      if (show.progress) cat("\t", "Max ICL : ", (round(maxICL, 2)), "\t", "Current ICL : ", (round(res1$ICL, 2)))
    }
    resRec[[i]] <- res1
  }
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t", "Trial : ", 1, " of ", 1)
  par <- res1Best$par[[which.max(res1Best$logLik)]]

  res2Best <- tryCatch(fit.fmmstil(x, K, param = par, show.progress = FALSE, control = control),
    error = function(e) res1Best,
    warning = function(w) res1Best
  )

  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (res2Best$ICL > maxICL) maxICL <- res2Best$ICL
  
  if (show.progress) cat("\t", "Max ICL : ", (round(maxICL, 2)), "\t", "Current ICL : ", (round(res2Best$ICL, 2)))

  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}



#' Automatic model based clustering via fmmstil.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param numTrial.fun a function of K that returns the number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{record}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil(as.matrix(log(RiverFlow)))
cluster.fmmstil <- function(x, numTrial.fun, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (missing(numTrial.fun)) numTrial.fun <- function(K) 2^K
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method

  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  maxICL <- -Inf
  while (TRUE) {
    res <- cluster.fmmstil.K(x, K, numTrial.fun(K), init.cluster.method, init.param.method, show.progress = show.progress, control = control)
    if (show.progress) cat("\n")
    
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    resRec[[K]]$unrestricted <- list(res$unrestricted)
  
    if (maxICL >= max(res$restricted$ICL, res$unrestricted$ICL)) {
      resBest$time <- difftime(Sys.time(), startTime, units = "secs")
      return(list(res = resBest, Record = resRec))
    } else if (res$restricted$ICL > res$unrestricted$ICL) {
      resBest <- res$restricted
      maxICL <- res$restricted$ICL
      K <- K + 1
    } else {
      resBest <- res$unrestricted
      maxICL <- res$unrestricted$ICL
      K <- K + 1
    }
  }
}



#' Clustering using mstil given a fixed number of clusters in parallel.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param K positive integer, number of cluster.
#' @param ncore number of cpu core to be used in parallel. By default 1.
#' @param numTrial a positive integer, number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{restricted}{a list containing details of the best fitted fmmstil.r.}
#' \item{unrestricted}{a list containing details of the best fitted fmmstil.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.L.parallel(as.matrix(log(RiverFlow)),2,2)
cluster.fmmstil.K.parallel <- function(x, K, ncore = 1, numTrial = 1, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgConR" %in% names(control)) control$cvgConR <- 1e-1

  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }

  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  
  if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Number of Trials : ", numTrial)
  
  seed <- sample(1000000000, numTrial)
  
  fit.fmmstil.r.seed <- function(seed){
    set.seed(seed)
    init.cluster <- init.cluster.method(x, K)
    res1 <- fit.fmmstil.r(x, K, init.cluster = init.cluster, init.param.method = init.param.method, show.progress = FALSE, control = control)
    par <- res1$par[[which.max(res1$logLik)]]
    fitness1 <- fmmstil.r.fitness(x, par)
    res1$ICL <- fitness1$ICL
    res1$clust <- fitness1$clust
    res1$AIC <- fitness1$AIC
    res1$BIC <- fitness1$BIC
    return(res1)
  }
  
  resRec <- parallel::mclapply(seed, fit.fmmstil.r.seed , mc.cores = ncore)

  
  ICLRec <- c()
  for (i in 1:numTrial) ICLRec <- c(ICLRec, resRec[[i]]$ICL)
  res1Best <- resRec[[which.max(ICLRec)]]
  par <- res1Best$par[[which.max(res1Best$logLik)]]

  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec), 2)))


  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
  res2Best <- tryCatch(fit.fmmstil(x, K, param = par, show.progress = FALSE, control = control),
    error = function(e) res1Best,
    warning = function(w) res1Best
  )

  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec, res2Best$ICL), 2)))
  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}



#' Automatic model based clustering via fmmstil in parallel.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param ncore number of cpu core to be used in parallel. By default 1.
#' @param numTrial.fun a function of K that returns the number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.parallel(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.parallel <- function(x, ncore = 1, numTrial.fun, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (missing(numTrial.fun)) numTrial.fun <- function(K) 2^K
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method

  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }


  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  maxICL <- -Inf
  while (TRUE) {
    res <- cluster.fmmstil.K.parallel(x, K, ncore, numTrial.fun(K), init.cluster.method, init.param.method, show.progress = show.progress, control = control)
    if (show.progress) cat("\n")
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    resRec[[K]]$unrestricted <- list(res$unrestricted)
    if (maxICL >= max(res$restricted$ICL, res$unrestricted$ICL)) {
      resBest$time <- difftime(Sys.time(), startTime, units = "secs")
      return(list(res = resBest, record = resRec))
    } else if (res$restricted$ICL > res$unrestricted$ICL) {
      resBest <- res$restricted
      maxICL <- res$restricted$ICL
      K <- K + 1
    } else {
      resBest <- res$unrestricted
      maxICL <- res$unrestricted$ICL
      K <- K + 1
    }
  }
}



#' This function returns the log likelihood, ICL, AIC, BIC, and group memberships for fmmstil
#' @param x a n x k matrix, representing n k-variate samples.
#' @param param list of lists of fitted parameters, containing omega, lambda, delta, Ainv, and nu.
#' @param u (Optional) a m x k matrix, a set of samples generated from standard k-dimensional multivariate t distribution with degree of freedom nu.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{number of samples used to estimate the density and log-likelihood functions. By default 1e6.}
##' }
#' @return a list with components:
#' \item{logLik}{the estimated log likelihood.}
#' \item{ICL}{the estiamted integrated completed log likelihood.}
#' \item{AIC}{the estiamted Akaike information criteria.}
#' \item{BIC}{the estiamted Bayesian information criteria.}
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



#' This function returns the log likelihood, ICL, AIC, BIC, and group memberships for fmmstil.r
#' @param x a n x k matrix, representing n k-variate samples.
#' @param param list of lists of fitted parameters, containing omega, lambda, delta, Ainv, and nu.
#' @return a list with components:
#' \item{logLik}{the log likelihood.}
#' \item{ICL}{the integrated completed log likelihood.}
#' \item{AIC}{the Akaike information criteria.}
#' \item{BIC}{the Bayesian information criteria.}
#' \item{clust}{the estimated group membership.}
#' @export
fmmstil.r.fitness <- function(x, param) {
  weight <- .fmmstil.r.weight(x, param$omega, param$lambda,param$delta, param$Ainv, param$nu)
  K <- length(param$omega)
  k <- ncol(x)
  m <- k * ( k + 1) / 2 + 1 + k + k
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



#' Clustering using mstil given a fixed number of clusters in parallel for windows. 
#' @param x a n x k matrix, representing n k-variate samples.
#' @param K positive integer, number of cluster.
#' @param ncore number of cpu core to be used in parallel. By default 1.
#' @param numTrial a positive integer, number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{restricted}{a list containing details of the best fitted fmmstil.r.}
#' \item{unrestricted}{a list containing details of the best fitted fmmstil.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.L.parallel(as.matrix(log(RiverFlow)),2,2)
cluster.fmmstil.K.parallel.windows <- function(x, K, ncore = 1, numTrial = 1, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgConR" %in% names(control)) control$cvgConR <- 1e-1
  if (ncore > parallel::detectCores()) {
    ncore <- parallel::detectCores()
    warning("Not enough available core")
  }
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Number of Trials : ", numTrial)
  
  ncore = min(ncore,numTrial)
  
  fit.fmmstil.r.seed.windows <- function(seed, data, init.cluster.method, init.param.method, control ){
    set.seed(seed)
    init.cluster = init.cluster.method(data,K)
    res1 <- fit.fmmstil.r(data, K, init.cluster = init.cluster, init.param.method = init.param.method, show.progress = FALSE, control = control)
    par <- res1$par[[which.max(res1$logLik)]]
    fitness1 <- fmmstil.r.fitness(x, par)
    res1$ICL <- fitness1$ICL
    res1$clust <- fitness1$clust
    res1$AIC <- fitness1$AIC
    res1$BIC <- fitness1$BIC
    return(res1)
  }
  
  seed <- sample(1000000000, numTrial)
  
  cl <- parallel::makePSOCKcluster(ncore)
  parallel::setDefaultCluster(cl)
  resRec <- parallel::parLapply(NULL, seed, fit.fmmstil.r.seed.windows, control = control, data = x, init.cluster.method = init.cluster.method, init.param.method = init.param.method)
  parallel::stopCluster(cl)
  ICLRec <- c()
  for (i in 1:numTrial) ICLRec <- c(ICLRec, resRec[[i]]$ICL)
  res1Best <- resRec[[which.max(ICLRec)]]
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec), 2)))
  
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
  res2Best <- tryCatch(fit.fmmstil(x, K, param = par, show.progress = FALSE, control = control),
                       error = function(e) res1Best,
                       warning = function(w) res1Best
  )
  
  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec, res2Best$ICL), 2)))
  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}


#' Automatic model based clustering via fmmstil in parallel for windows.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param ncore number of cpu core to be used in parallel. By default 1.
#' @param numTrial.fun a function of K that returns the number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgConR is 0.1 instead.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.parallel(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.parallel.windows <- function(x, ncore = 1, numTrial.fun, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (missing(numTrial.fun)) numTrial.fun <- function(K) 2^K
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  
  
  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  maxICL <- -Inf
  while (TRUE) {
    res <- cluster.fmmstil.K.parallel.windows(x, K, ncore, numTrial.fun(K), init.cluster.method, init.param.method, show.progress = show.progress, control = control)
    if (show.progress) cat("\n")
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    resRec[[K]]$unrestricted <- list(res$unrestricted)
    if (maxICL >= max(res$restricted$ICL, res$unrestricted$ICL)) {
      resBest$time <- difftime(Sys.time(), startTime, units = "secs")
      return(list(res = resBest, record = resRec))
    } else if (res$restricted$ICL > res$unrestricted$ICL) {
      resBest <- res$restricted
      maxICL <- res$restricted$ICL
      K <- K + 1
    } else {
      resBest <- res$unrestricted
      maxICL <- res$unrestricted$ICL
      K <- K + 1
    }
  }
}



