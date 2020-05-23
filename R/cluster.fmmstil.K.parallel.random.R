#' Clustering using mstil given a fixed number of clusters in parallel using random initial parameters. 
#' @param x matrix of quantiles of size n x k. Each row is taken as a quantile.
#' @param K a positive integer, represents the number of clusters in the model. 
#' @param ncore a positive integer, represents the number of cpu threads to be used in parallel. By default 1.
#' @param numTrial a positive integer, represents the number of trials to be evaluated. By default 1.
#' @param criteria Either 'ICL', 'BIC', or 'AIC'. Represents the type of information criteria used for model selection. By default 'ICL'.
#' @param init.cluster.method a function of x, K that seperates x into K initial clusters.
#' @param init.param.method a function of x, returns initial parameters.
#' @param show.progress a logical value. If TRUE, progress of the algorithm will be printed in console. By default TRUE.
#' @param control list of control variables, it accepts all control arguments used in fit.fmmstil.r and fit.fmmsil. In this case, the default cvgTolR is 1e-1.
#' @return a list with components:
#' \item{restricted}{a list containing details of the best fitted fmmstil.r.}
#' \item{unrestricted}{a list containing details of the best fitted fmmstil.}
#' \item{recordR}{a list of list containing details all fitted fmmstil.r.}
#' @export
#' @examples
#' # Not run:
#' # data(RiverFlow)
#' # cluster.fmmstil.K.parallel.random(as.matrix(log(RiverFlow)),2,2)
cluster.fmmstil.K.parallel.random <- function(x, K, ncore = 1, numTrial = 1, criteria = c('ICL', 'BIC', 'AIC'), init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  .check.control(control)
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-2
  if (ncore > parallel::detectCores()) {
    ncore <- parallel::detectCores()
    warning("Not enough available core")
  }

  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method.random
  if (missing(init.param.method)) init.param.method <- .default.init.param.method.random
  if (missing(criteria)) criteria <- 'ICL'
  if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Number of Trials : ", numTrial)
  seedFmmstil <- sample(.Machine$integer.max, 1)
  ncore1 = min(ncore,numTrial)
  
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
  
  seed <- sample(.Machine$integer.max, numTrial)
  fit.fmmstil.r.seed.windows <- function(trial, data, K, initParamList, seed, control = list()){
    set.seed(seed[trial])
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
  
  
  cl <- parallel::makePSOCKcluster(ncore1)
  parallel::setDefaultCluster(cl)
  resRec <- parallel::parLapply(NULL, 1:numTrial, fit.fmmstil.r.seed.windows, seed = seed, control = control, data = x, K = K, initParamList = initParamList)
  parallel::stopCluster(cl)
  
  set.seed(seedFmmstil)
  
  criteriaRec <- c()
  for (i in 1:numTrial) criteriaRec <- c(criteriaRec, resRec[[i]][[criteria]])
  res1Best <- resRec[[which.min(criteriaRec)]]
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  
  if (show.progress) cat("\t", "Min. ", criteria, " : ", (round(min(criteriaRec), 2)))
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
  res2Best <- tryCatch(fit.fmmstil.parallel(x, K, ncore = ncore, param = par, show.progress = FALSE, control = control),
                       error = function(e) return(res1Best))
  
  par <- res2Best$par[[which.max(res2Best$logLik)]]
  fitness2 <- fmmstil.fitness(x, par, control = control)
  res2Best$ICL <- fitness2$ICL
  res2Best$clust <- fitness2$clust
  res2Best$AIC <- fitness2$AIC
  res2Best$BIC <- fitness2$BIC
  
  if (show.progress) cat("\t", "Min. ", criteria, " : ", (round(min(criteriaRec, res2Best[[criteria]]), 2)))
  
  
  return(list(restricted = res1Best, unrestricted = res2Best, recordR = resRec))
}
