#' Automatic model based clustering via fmmstil in parallel using divisive hierarchical method.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param ncore a positive integer, represents the number of cpu threads to be used in parallel. By default 1.
#' @param criteria Either 'ICL', 'BIC', or 'AIC'. Represents the type of information criteria used for model selection. By default 'ICL'.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a function of x, returns initial parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 1e-2 and the default cvgTolR is 1e-1.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{record}{a list of lists containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.parallel.divisive(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.parallel.divisive <- function(x, ncore = 1, criteria = c('ICL', 'BIC', 'AIC'), init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (missing(init.cluster.method))init.cluster.method <- .default.init.cluster.method.kmeans
  if (missing(init.param.method)) init.param.method <- .default.init.param.method.t
  if (missing(criteria)) criteria <- 'ICL'
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  cluster0 <- rep(1, nrow(x))
  minCritiera <- Inf
  while (TRUE) {
    if (K > 2) cluster0 <- res$unrestricted$clust
    if (K == 1) res <- cluster.fmmstil.K.parallel.random(x, K, ncore, criteria = criteria, init.cluster.method = init.cluster.method, init.param.method = init.param.method, control = control, show.progress = show.progress)
    else res <- .cluster.fmmstil.K.parallel.divisive(x = x, K = K, ncore = ncore, cluster0 = cluster0, criteria = criteria, init.cluster.method = init.cluster.method, init.param.method = init.param.method, show.progress = show.progress, control = control)
    if (show.progress) cat("\n")
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    resRec[[K]]$unrestricted <- list(res$unrestricted)
    if (minCritiera <= min(res$restricted[[criteria]], res$unrestricted[[criteria]])) {
      resBest$time <- difftime(Sys.time(), startTime, units = "secs")
      return(list(res = resBest, record = resRec))
    } else if (res$restricted[[criteria]] < res$unrestricted[[criteria]]) {
      resBest <- res$restricted
      minCritiera <- res$restricted[[criteria]]
      K <- K + 1
    } else {
      resBest <- res$unrestricted
      minCritiera <- res$unrestricted[[criteria]]
      K <- K + 1
    }
  }
}
