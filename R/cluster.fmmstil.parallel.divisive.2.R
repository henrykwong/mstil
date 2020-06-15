#' Automatic model based clustering via fmmstil in parallel using divisive hierarchical method 2. 
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
#' # cluster.fmmstil.parallel.divisive.2(as.matrix(log(RiverFlow)))
cluster.fmmstil.parallel.divisive.2 <- function(x, ncore = 1, criteria = c('ICL', 'BIC', 'AIC'), init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (missing(init.cluster.method))init.cluster.method <- .default.init.cluster.method.kmeans
  if (missing(init.param.method)) init.param.method <- .default.init.param.method.t
  if (missing(criteria)) criteria <- 'ICL'
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  startTime <- Sys.time()
  
  
  res1 <- cluster.fmmstil.r.parallel.divisive(x, ncore, criteria, 
                                              init.cluster.method, init.param.method,
                                              show.progress, control)
  
  res1Best <- res1$res
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  K = length(par$omega)
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
  res2Best <- tryCatch(fit.fmmstil.parallel(x, K - 1, ncore = ncore, param = par, show.progress = FALSE, control = control)
                       , error = function(e) return(res1Best))
  
  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par, control = control)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (show.progress) cat("\t", "Min. ", criteria, " : ", (round(min(res1Best[[criteria]], res2Best[[criteria]]), 2)))
  
  if (res2Best[[criteria]] < res1Best[[criteria]]){
    resBest <- res2Best
  } else resBest <- res1Best

  return(list(res = resBest, recordR = res1$recordR))
}
