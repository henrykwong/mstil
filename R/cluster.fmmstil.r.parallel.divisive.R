#' Automatic model based clustering via fmmstil.r in parallel using divisive hierarchical method.
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param ncore a positive integer, represents the number of cpu threads to be used in parallel. By default 1.
#' @param criteria Either 'ICL', 'BIC', or 'AIC'. Represents the type of information criteria used for model selection. By default 'ICL'.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a function of x, returns initial parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil.
#' @return a list with components:
#' \item{res}{a list containing details of the best fitted distribution.}
#' \item{recordR}{a list of lists containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.r.parallel.divisive(as.matrix(log(RiverFlow)))
cluster.fmmstil.r.parallel.divisive <- function(x, ncore = 1, criteria = c('ICL', 'BIC', 'AIC'), init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
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
  minCriteria <- Inf
  
  while (TRUE) {
    
    if (K > 2) cluster0 <- res$res$clust
    if (K == 1){
      if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", 1, "\t", "Number of Trials : ", 1)
      res <- fit.fmmstil.r(x, 1, show.progress = FALSE, control = control)
      par <- res$par[[which.max(res$logLik)]]
      fitness1 <- fmmstil.r.fitness(x, par)
      res$ICL <- fitness1$ICL
      res$clust <- fitness1$clust
      res$AIC <- fitness1$AIC
      res$BIC <- fitness1$BIC
      if (show.progress) cat("\t", "Min. ", criteria, " : ", (round(res[[criteria]], 2)))
      res = list(res = res, recordR = list(res))
    }
    else{
      res <- .cluster.fmmstil.r.K.parallel.divisive(x = x, K = K, ncore = ncore, 
                                                       cluster0 = cluster0, criteria = criteria, 
                                                       init.cluster.method = init.cluster.method, 
                                                       init.param.method = init.param.method, 
                                                       show.progress = show.progress, control = control)
    }
    
    resRec[[K]] <- list()
    resRec[[K]]$restricted <- res$recordR
    if (minCriteria <= res$res[[criteria]]) {
      break
    } else {
      resBest <- res$res
      minCriteria <- res$res[[criteria]]
      K <- K + 1
    }
  }
  return(list(res = resBest, recordR = resRec))
}
