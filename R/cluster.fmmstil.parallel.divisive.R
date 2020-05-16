#' Automatic model based clustering via fmmstil in parallel using divisive hierarchical method.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param ncore a positive integer, represents the number of cpu threads to be used in parallel. By default 1.
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
#' # cluster.fmmstil.parallel.divisive(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.parallel.divisive <- function(x, ncore = 1, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (missing(init.cluster.method))init.cluster.method <- .default.init.cluster.method2
  if (missing(init.param.method)) init.param.method <- .default.init.param.method2
  
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  
  resRec <- list()
  startTime <- Sys.time()
  K <- 1
  cluster0 <- rep(1, nrow(x))
  maxICL <- -Inf
  while (TRUE) {
    if (K > 2) cluster0 <- res$unrestricted$clust
    if (K == 1) res <- cluster.fmmstil.K.parallel.random(x, K, ncore, init.cluster.method = init.cluster.method, init.param.method = init.param.method, control = control, show.progress = show.progress)
    else res <- .cluster.fmmstil.K.parallel.divisive(x, K, ncore, cluster0 = cluster0, init.cluster.method, init.param.method, show.progress = show.progress, control = control)
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
