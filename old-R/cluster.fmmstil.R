#' Automatic model based clustering via fmmstil.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param numTrial.fun a function of K that returns the number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgTolR is 0.1 instead.
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
