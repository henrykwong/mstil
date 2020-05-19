#' Automatic model based clustering via fmmstil in parallel using random initial parameters
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param ncore a positive integer, represents the number of cpu threads to be used in parallel. By default 1.
#' @param numTrial.fun a function of K that returns the number of trials to be evaluated.
#' @param criteria Either 'ICL', 'BIC', or 'AIC'. Represents the type of information criteria used for model selection. By default 'ICL'.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters. 
#' @param init.param.method a function of x, returns initial parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 1e-2 and the default cvgTolR is 1e-1.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{recordR}{a list of lists containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.parallel.random(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.parallel.random <- function(x, ncore = 1, numTrial.fun, criteria = c('ICL', 'BIC', 'AIC'), init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (missing(numTrial.fun)) numTrial.fun <- function(K) 2^K
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method.random
  if (missing(init.param.method)) init.param.method <- .default.init.param.method.random
  if (missing(criteria)) criteria <- 'ICL'
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  
  
  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  minCriteria <- Inf
  while (TRUE) {
    res <- cluster.fmmstil.K.parallel.random(x, K, ncore, numTrial.fun(K), init.cluster.method, init.param.method, show.progress = show.progress, control = control)
    if (show.progress) cat("\n")
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    resRec[[K]]$unrestricted <- list(res$unrestricted)
    if (minCriteria <= min(res$restricted[[criteria]], res$unrestricted[[criteria]])) {
      resBest$time <- difftime(Sys.time(), startTime, units = "secs")
      return(list(res = resBest, record = resRec))
    } else if (res$restricted[[criteria]] < res$unrestricted[[criteria]]) {
      resBest <- res$restricted
      minCriteria <- res$restricted[[criteria]]
      K <- K + 1
    } else {
      resBest <- res$unrestricted
      minCriteria <- res$unrestricted[[criteria]]
      K <- K + 1
    }
  }
}
