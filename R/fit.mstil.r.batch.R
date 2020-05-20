#' This function finds the maximum penalised quasi likelihood estiamtes for mstil.r using batch of subsamples.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param param list of inital parameters, contains lambda, delta, Ainv, and nu.
#' @param show.progress logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, see 'details'.
#' @details The control argument is a list that accepts the following components.
##' \describe{
##'  \item{cvgNR}{a positive integer. The algorithm stops when the estimated log-likelihood is not improved by at least cvgTolR in cvgNR iterations. By default 5.}
##'  \item{cvgTolR}{a positive value. The algorithm stops when the estimated log-likelihood is not improved by at least cvgTolR in cvgNR iterations. By default 1e-2.}
##'  \item{lambdaPenalty}{a positive value, represents the L1 penalty coefficient for lambda. By default 0.}
##'  \item{ainvPenalty}{a positive value, represents the L2 penalty coefficient for Ainv. By default 0.}
##'  \item{maxitR}{a positive integer, represents the maximum number iterations allowed. By default 1e3.}
##'  \item{maxitOptimR}{a positive integer, represents the maximum number of iterations allowed in optim. By default 10.}
##'  \item{batchSizeR}{a positive integer, represents the batch sample size. By default n.}
##' }
#' @return a list with components:
#' \item{logLik}{a vector of values of the log-likelihood function after each itereation.}
#' \item{par}{a list of lists of fitted parameters after each iteration.}
#' \item{time}{a vector recorded the time elapsed after each iteration.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # fit.mstil.r.batch(as.matrix(log(RiverFlow)), control = list(batchSizeR = 100))
fit.mstil.r.batch <- function(x, param, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (!"maxitR" %in% names(control)) control$maxitR <- 1e3
  if (!"cvgNR" %in% names(control)) control$cvgNR <- 5
  if (!"batchSizeR" %in% names(control)) control$batchSizeR <- nrow(x)
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-2
  maxitR <- control$maxitR
  cvgNR <- control$cvgNR
  cvgTolR <- control$cvgTolR
  
  batchSizeR <- min(nrow(x), control$batchSizeR)

  if (missing(param)) param <- .default.init.param.method.t(x)
  .check.mstil.r.param(ncol(x), param$lambda, param$delta, param$Ainv, param$nu)
  
  res <- list(lambda = param$lambda, delta = param$delta, Ainv = param$Ainv, nu = param$nu)
  
  startTime <- Sys.time()
  
  lik <- sum(dmstil.r(x, res$lambda, res$delta, res$Ainv, res$nu, log.p = TRUE))
  likRec <- lik
  timeRec <- 0
  resRec <- list()
  resRec[[1]] <- res
  
  
  for (i in 2:(maxitR + 1)) {
    xBatch <- x[sample(nrow(x), batchSizeR),]
    res <- fit.mstil.r(xBatch, res, control = control)
    lik <- sum(dmstil.r(x, res$lambda, res$delta, res$Ainv, res$nu, log.p = TRUE))
    resRec[[i]] <- res
    likRec <- c(likRec, lik)
    timeRec <- c(timeRec, as.numeric(Sys.time() - startTime, units = "secs"))
    if (show.progress) {
      cat("\r", "Iteration : ", (i - 1), "Current Likelihood : ", round(likRec[length(likRec)]), "Maximum Likelihood : ", round(max(likRec)), "\t")
    }
    
    if (i > cvgNR) {
      if (all((likRec[i - cvgNR] + cvgTolR) > likRec[(i - cvgNR + 1):i])) {
        if (show.progress) cat("\n", "Converged!")
        return(list(logLik = likRec, par = resRec, time = timeRec))
        break
      }
    }
  }
  if (show.progress) cat("\n", "Maximum number of iterations reached!")
  return(list(logLik = likRec, par = resRec, time = timeRec))
}
