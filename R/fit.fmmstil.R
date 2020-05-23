#' This function finds the penalised quasi likelihood estiamtes for finite mixture of mstil via EM.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param K number of clusters.
#' @param param list of lists of inital parameters, contains omega, lambda, delta, Ainv, and nu.
#' @param init.cluster (optional) initial clusters used to find initial parameters.
#' @param init.param.method (optional) method to obtain initial parameters. It needs to be a function of x and K that returns a list of lists of appropriate parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{a positive integer, represents the number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{a value between 0.5 and 1, represents the confidence level of the log-likelihood to be calculated. By default 0.95.}
##'  \item{cvgN}{a positive integer. The algorithm stops when the estimated log-likelihood is not improved in cvgN iterations. By default 5.}
##'  \item{lambdaPenalty}{a positive value, represents the L2 penalty coefficient for lambda. By default 1e-4.}
##'  \item{ainvPenalty}{a positive value, represents the L2 penalty coefficient for Ainv. By default 1e-6.}
##'  \item{maxit}{a positive integer, represents the maximum number of EM iterations allowed. By default 1e3.}
##'  \item{maxitOptim}{a positive integer, represents the maximum number of iterations in optim allowed within each M-step. By default 10.}
##'  \item{numGradSample}{a positive integer, represents the number of samples used to estimate the gradient. By default 1e4.}
##'  \item{finDiffStep}{a positive value, represents the step size to estimate gradient w.r.t. nu. By default 1e-5.}
##'  \item{stepSgd}{a positive value, represents the step size to be used in the stochastic gradient step to optimise nu. By default 1e-2.}
##'  \item{iterSgd}{a positive integer, represents the number of iterations to be used in the stochastic gradient step to optimise nu. By default 1e2.}
##'  \item{dimRateSgd}{a positive value, represents the diminishing rate for step size to be used in the stochastic gradient step to optimise nu. By default 1e-2.}
##'  \item{batchSize}{a positive integer, represents the batch sample size. By default n.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector of the estimated log-likelihood after each itereation.}
#' \item{par}{a list of lists of lists of fitted parameters after each iteration.}
#' \item{time}{a vector recorded the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.fmmstil(as.matrix(log(RiverFlow)), 2)
fit.fmmstil <- function(x, K, param, init.cluster, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  
  if (!"maxit" %in% names(control)) control$maxit <- 1e3
  if (!"cvgN" %in% names(control)) control$cvgN <- 5
  if (!"batchSize" %in% names(control)) control$batchSize <- nrow(x)
  maxit <- control$maxit
  cvgN <- control$cvgN
  batchSize <- min(control$batchSize, nrow(x))
  
  n <- nrow(x)
  
  if (missing(param)) {
    param <- list(omega = list(), lambda = list(), delta = list(), Ainv = list(), nu = list())
    if (missing(init.param.method)) init.param.method <- .default.init.param.method.random
    if (missing(init.cluster)) init.cluster <- .default.init.cluster.method.random(x, K)
    param$omega <- as.list(table(init.cluster) / n)
    for (i in 1:K) {
      initFit <- init.param.method(x[which(init.cluster == unique(init.cluster)[i]), ])
      param$lambda[[i]] <- initFit$lambda
      param$delta[[i]] <- initFit$delta
      param$Ainv[[i]] <- initFit$Ainv
      param$nu[[i]] <- initFit$nu
    }
  }
  .check.fmmstil.param(ncol(x), param)
  
  res <- list(omega = param$omega, lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  res1 <- res
  startTime <- Sys.time()
  likRec <- c()
  resRec <- list()
  timeRec <- c()
  for (i in 1:maxit) {
    w_1 <- .fmmstil.weight(x, res1$omega, res1$lambda, res1$delta, res1$Ainv, res1$nu, control = control)
    d <- rowSums(w_1)
    logLik <- sum(log(d))
    if (!is.na(logLik)){
      w_ <- w_1
      w <- w_ / d
      res <- res1
    }
    res$omega <- as.list(colSums(w) / n)
    likRec <- c(likRec, logLik)
    resRec[[i]] <- res
    timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))
    if (i > cvgN) {
      if (all(likRec[i - cvgN] > likRec[(i - cvgN + 1):i])) {
        if (show.progress) cat("\n", "converged!")
        return(list(logLik = likRec, par = resRec, time = timeRec))
        break
      }
    }
    if (show.progress) {
      cat("\r", "Iteration : ", i, "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    batch <- sample(nrow(x),batchSize)
    for (j in 1:K) {
      res2 <- .fit.mstil.1.weighted(x[batch,], w[batch, j], lambda = res$lambda[[j]], delta = res$delta[[j]], Ainv = res$Ainv[[j]], nu = res$nu[[j]], control = control)
      res1$lambda[[j]] <- res2$lambda
      res1$delta[[j]] <- res2$delta
      res1$Ainv[[j]] <- res2$Ainv
      res1$nu[[j]] <- .fit.mstil.2.weighted(x[batch,], w[batch, j], lambda = res1$lambda[[j]], delta = res1$delta[[j]], Ainv = res1$Ainv[[j]], nu = res1$nu[[j]], control = control)
    }
  }
  
  w_1 <- .fmmstil.weight(x, res1$omega, res1$lambda, res1$delta, res1$Ainv, res1$nu, control = control)
  d <- rowSums(w_1)
  logLik <- sum(log(d))
  if (!is.na(logLik)){
    w_ <- w_1
    w <- w_ / d
    res <- res1
  } else{
    logLik <- likRec[length(likRec)]
  }
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
