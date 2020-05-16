#' Clustering using mstil given a fixed number of clusters in parallel.
#' @param x a n x k matrix, representing n k-variate samples.
#' @param K positive integer, number of cluster.
#' @param ncore number of cpu core to be used in parallel. By default 1.
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
#' # cluster.fmmstil.L.parallel(as.matrix(log(RiverFlow)),2,2)
cluster.fmmstil.K.parallel <- function(x, K, ncore = 1, numTrial = 1, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-1
  
  if (ncore > parallel::detectCores()) {
    warning("Not enough available core")
    ncore <- parallel::detectCores()
  }
  
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method
  if (missing(init.param.method)) init.param.method <- .default.init.param.method
  
  if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Number of Trials : ", numTrial)
  
  seedFmmstil <- sample(.Machine$integer.max, 1)
  
  initParamList <- list()
  for (trial in 1:numTrial){
    initParamList[[trial]] <- list(omega = list(), lambda = list(), delta = list(), Ainv = list(), nu = list())
    initCluster <- init.cluster.method(x, K)
    initParamList[[trial]]$omega <- as.list(table(initCluster) / length(initCluster))
    for (kk in 1:K){
      initFit <- init.param.method(x[which(initCluster == kk),])
      initParamList[[trial]]$lambda[[kk]] = initFit$lambda
      initParamList[[trial]]$delta[[kk]] = initFit$delta
      initParamList[[trial]]$Ainv[[kk]] = initFit$Ainv
      initParamList[[trial]]$nu[[kk]] = initFit$nu
    }
  }
  
  
  fit.fmmstil.r.seed <- function(trial, data, K, initParamList, control = list()){
    init.cluster = init.cluster.method(data,K)
    res1 <- fit.fmmstil.r(data, K, initParamList[[trial]], show.progress = FALSE, control = control)
    par <- res1$par[[which.max(res1$logLik)]]
    fitness1 <- fmmstil.r.fitness(x, par)
    res1$ICL <- fitness1$ICL
    res1$clust <- fitness1$clust
    res1$AIC <- fitness1$AIC
    res1$BIC <- fitness1$BIC
    return(res1)
  }
  
  
  resRec <- parallel::mclapply(1:numTrial, fit.fmmstil.r.seed , mc.cores = ncore, data = x, K = K, initParamList = initParamList, control = control)
  
  set.seed(seedFmmstil)
  ICLRec <- c()
  for (i in 1:numTrial) ICLRec <- c(ICLRec, resRec[[i]]$ICL)
  res1Best <- resRec[[which.max(ICLRec)]]
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec), 2)))
  
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
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
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec, res2Best$ICL), 2)))
  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}
