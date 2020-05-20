#' This function finds the maximum penalised quasi likelihood estiamtes for mstil.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of inital parameters, contains lambda, delta, Ainv, and nu.
#' @param show.progress logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{numLikSample}{a positive integer, represents the number of samples used to estimate the density and log-likelihood functions. By default 1e6. }
##'  \item{conLevel}{a value between 0.5 and 1, represents the confidence level of the log-likelihood to be calculated. By default 0.95.}
##'  \item{cvgN}{a positive interger. The algorithm stops when the estimated log-likelihood is not improved in cvgN iterations. By default 5.}
##'  \item{lambdaPenalty}{a positive value, represents the L1 penalty coefficient for lambda. By default 0.}
##'  \item{ainvPenalty}{a positive value, represents the L2 penalty coefficient for Ainv. By default 0.}
##'  \item{maxit}{a positive integer, represents the maximum number iterations allowed. By default 0.}
##'  \item{maxitOptim}{maximum number of iterations in optim. By default 10.}
##'  \item{numGradSample}{a positive integer, represents the number of samples used to estimate the gradient. By default 1e4.}
##'  \item{finDiffStep}{a positive value, represents the step size used to estimate gradient w.r.t. nu. By default 1e-5.}
##'  \item{stepSgd}{a positive value, represents the step size to be used in the stochastic gradient step to optimise nu. By default 1e-2.}
##'  \item{iterSgd}{a positive integer, represents the number of iterations to be used in the stochastic gradient step to optimise nu. By default 1e2.}
##'  \item{dimRateSgd}{a positive value, represents the diminishing rate for step size to be used in the stochastic gradient step to optimise nu. By default 1e-2.}
##'  \item{batchSize}{a positive integer, represents the batch sample size. By default n.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector of estimated values of the log-likelihood function after each itereation.}
#' \item{par}{a list of lists of fitted parameters after each iteration.}
#' \item{logLikLower}{a vector of the lower bound of the estimated log-likelihood function after each iteration.}
#' \item{logLikUpper}{a vector of the upper bound of the estimated log-likelihood function after each iteration.}
#' \item{time}{a vector recorded the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.mstil(as.matrix(log(RiverFlow)))
fit.mstil <- function(x, param, show.progress = TRUE, control = list()) {
  .check.control(control)
  
  if (!"maxit" %in% names(control)) control$maxit <- 1e3
  if (!"cvgN" %in% names(control)) control$cvgN <- 5
  if (!"batchSize" %in% names(control)) control$batchSize <- nrow(x)
  maxit <- control$maxit
  cvgN <- control$cvgN
  batchSize <- min(control$batchSize, nrow(x))
  
  if (missing(param)){
    if (batchSize == nrow(x)) param <- fit.mstil.r(x, control = control)
    else{
      control$batchSizeR <- control$batchSize
      initFit <- fit.mstil.r.batch(x, show.progress = FALSE, control = control)
      param <- initFit$par[[which.max(initFit$logLik)]]
    }
  }
  .check.mstil.param(ncol(x), param$lambda, param$delta, param$Ainv, param$nu)
  
  res <- list(lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  startTime <- Sys.time()
  lik <- mstil.logLik(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
  likRec <- lik$logLik
  likLowerRec <- lik$logLikLower
  timeRec <- 0
  resRec <- list()
  resRec[[1]] <- res
  
  for (i in 2:(maxit + 1)) {
    xBatch <- x[sample(nrow(x), batchSize),]
    res <- .fit.mstil.1(xBatch, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    res$nu <- .fit.mstil.2(xBatch, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    lik <- mstil.logLik(x, res$lambda, res$delta, res$Ainv, res$nu, control = control)
    resRec[[i]] <- res
    likRec <- c(likRec, lik$logLik)
    likLowerRec <- c(likLowerRec, lik$logLikLowerBound)
    timeRec <- c(timeRec, as.numeric(Sys.time() - startTime, units = "secs"))
    if (show.progress) {
      cat("\r", "Iteration : ", (i - 1), "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    if (i > cvgN) {
      if (all(likRec[i - cvgN] > likRec[(i - cvgN + 1):i])) {
        if (show.progress) cat("\n", "Converged!")
        return(list(logLik = likRec, par = resRec, logLikLower = likLowerRec, time = timeRec))
        break
      }
    }
  }
  if (show.progress) cat("\n", "Maximum number of iterations reached!")
  return(list(logLik = likRec, par = resRec, logLikLowerBound = likLowerRec, time = timeRec))
}
