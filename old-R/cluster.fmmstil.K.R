#' Clustering using mstil given a fixed number of clusters.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param K positive integer, number of cluster.
#' @param numTrial a positive integer, number of trials to be evaluated.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a functino of x, return initial parameters.
#' @param show.progress show progress on console.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default lambdaPenalty is 0.01 and the default cvgTolR is 0.1 instead.
#' @return a list with components:
#' \item{restricted}{a list containing details of the best fitted fmmstil.r.}
#' \item{unrestricted}{a list containing details of the best fitted fmmstil.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.K(as.matrix(log(RiverFlow)),2)
cluster.fmmstil.K <- function(x, K, numTrial, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-1
  
  
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  
  resRec <- list()
  
  maxICL <- -Inf
  for (i in 1:numTrial) {
    if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Trial : ", i, " of ", numTrial)
    
    init.cluster <- init.cluster.method(x, K)
    
    res1 <- tryCatch(fit.fmmstil.r(x, K, init.cluster = init.cluster, init.param.method = init.param.method, show.progress = FALSE, control = control),
                     error = function(e) NA,
                     warning = function(w) NA
    )
    if (is.list(res1)) {
      par <- res1$par[[which.max(res1$logLik)]]
      fitness1 <- fmmstil.r.fitness(x, par)
      res1$ICL <- fitness1$ICL
      res1$clust <- fitness1$clust
      res1$AIC <- fitness1$AIC
      res1$BIC <- fitness1$BIC
      
      if (res1$ICL > maxICL) {
        res1Best <- res1
        maxICL <- res1$ICL
        guess1Best <- res1$clust
      }
      
      if (show.progress) cat("\t", "Max ICL : ", (round(maxICL, 2)), "\t", "Current ICL : ", (round(res1$ICL, 2)))
    }
    resRec[[i]] <- res1
  }
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t", "Trial : ", 1, " of ", 1)
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  
  res2Best <- tryCatch(fit.fmmstil(x, K, param = par, show.progress = FALSE, control = control),
                       error = function(e) res1Best,
                       warning = function(w) res1Best
  )
  
  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (res2Best$ICL > maxICL) maxICL <- res2Best$ICL
  
  if (show.progress) cat("\t", "Max ICL : ", (round(maxICL, 2)), "\t", "Current ICL : ", (round(res2Best$ICL, 2)))
  
  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}

