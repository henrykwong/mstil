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
##'  \item{cvgNR}{a positive integer. The algorithm stops when the estimated log-likelihood is not improved by at least cvgTolR in cvgNR iterations. By default 5.}
##'  \item{cvgTolR}{a positive value. The algorithm stops when the estimated log-likelihood is not improved by at least cvgTolR in cvgNR iterations. By default 1e-2.}
##'  \item{lambdaPenalty}{a positive value, represents the L2 penalty coefficient for lambda. By default 1e-4.}
##'  \item{ainvPenalty}{a positive value, represents the L2 penalty coefficient for Ainv. By default 1e-6.}
##'  \item{maxitR}{a positive integer, represents the maximum number of EM iterations allowed. By default 1000.}
##'  \item{maxitOptimR}{a positive integer, represents the maximum number of iterations in optim allowed within each M-step. By default 1e2.}
##'  \item{batchSizeR}{a positive integer, represents the batch sample size. By default n.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector of the estimated log-likelihood after each itereation.}
#' \item{par}{a list of lists of lists of fitted parameters after each iteration.}
#' \item{time}{a vector recorded the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.fmmstil.r(as.matrix(log(RiverFlow)), 2)
fit.fmmstil.r <- function(x, K, param, init.cluster, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (!"maxitR" %in% names(control)) control$maxitR <- 1e3
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-2
  if (!"cvgNR" %in% names(control)) control$cvgNR <- 5
  if (!"batchSizeR" %in% names(control)) control$batchSizeR <- nrow(x)
  maxitR <- control$maxitR
  cvgNR <- control$cvgNR
  cvgTolR <- cvgNR * control$cvgTolR

  n <- nrow(x)
  
  batchSizeR <- min(nrow(x), control$batchSizeR)
  
  
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
  
  .check.fmmstil.r.param(ncol(x), param)
  
  res <- list(omega = param$omega, lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  res1 <- res
  startTime <- Sys.time()
  likRec <- c()
  resRec <- list()
  timeRec <- c()
  for (i in 1:maxitR) {
    w_1 <- .fmmstil.r.weight(x, res1$omega, res1$lambda, res1$delta, res1$Ainv, res1$nu)
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
    resRec[[i]] <- res
    timeRec <- c(timeRec, difftime(Sys.time(), startTime, units = "secs"))
    
    
    if (i > cvgNR) {
      if (all((likRec[i - cvgNR] + cvgTolR) > likRec[(i - cvgNR + 1):i])) {
        if (show.progress) cat("\n", "Converged!")
        return(list(logLik = likRec, par = resRec, time = timeRec))
        break
      }
    }
    
    
    
    if (show.progress) {
      cat("\r", "Iteration : ", i, "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    for (j in 1:K) {
      batch <- sample(nrow(x),batchSizeR)
      res2 <- .fit.fmmstil.r.weighted(x[batch,], w[batch, j], lambda = res$lambda[[j]], delta = res$delta[[j]], Ainv = res$Ainv[[j]], nu = res$nu[[j]], control = control)
      res1$lambda[[j]] <- res2$lambda
      res1$delta[[j]] <- res2$delta
      res1$Ainv[[j]] <- res2$Ainv
      res1$nu[[j]] <- res2$nu
    }
  }
  
  w_1 <- .fmmstil.r.weight(x, res1$omega, res1$lambda, res1$delta, res1$Ainv, res1$nu)
  d <- rowSums(w_1)
  logLik <- sum(log(d))
  if (!is.na(logLik)){
    w_ <- w_1
    w <- w_ / d
    res <- res1
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

