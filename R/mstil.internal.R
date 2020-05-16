#' @keywords internal
.dmvt2 <- function(x, delta, Ainv, nu = 1, log.p = TRUE) {
  k <- ncol(x)
  xA <- Ainv %*% (t(x) - delta)
  xSx <- colSums(xA^2)
  if (nu > 1e6) {
    logDensity <- sum(log(abs(diag(Ainv)))) - 0.5 * k * log(2 * pi) - 0.5 * xSx
  } else {
    logDensity <- sum(log(abs(diag(Ainv)))) + lgamma((k + nu) / 2) - lgamma(nu / 2) -
      k / 2 * log(pi * nu) - 0.5 * (nu + k) * log1p(xSx / nu)
  }
  if (log.p) {
    return(logDensity)
  } else {
    return(exp(logDensity))
  }
}

#' @keywords internal
.mstil.grad.1 <- function(x, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  numGradSample <- control$numGradSample

  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  if (missing(u)) u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)
  m <- nrow(u)
  y <- t(t(x) - delta)
  xA <- Ainv %*% t(y)
  z <- t(xA)
  xSx <- colSums(xA^2)
  C <- 1 / (1 + 1 / nu * xSx)

  lambdaz <- z %*% lambda
  lambdau <- u %*% lambda
  Gz <- stats::plogis(lambdaz, log.p = TRUE)
  Gu <- stats::plogis(lambdau, log.p = TRUE)
  gz <- stats::dlogis(lambdaz, log = TRUE)
  gu <- stats::dlogis(lambdau, log = TRUE)

  gGz <- exp(gz - Gz)
  gGu <- exp(gu - Gu)
  expGu <- exp(rowSums(Gu))
  EGu <- mean(expGu)

  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(2 / nu * (y[, t2] * t(xA)[, t1]) * C)
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + n * diag(1 / diag(Ainv))
  dAinv1 <- grad
  dAinv2 <- t(t(matrix(colSums(y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)
  dmu1 <- rowSums(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C)))
  dmu2 <- as.vector(-(t(colSums(gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda1 <- matrix(colSums(gGz[, rep(1:p, each = k)] * t(xA[rep(1:k, p), ])), nrow = k)
  dlambda2 <- matrix(colSums(expGu * gGu[, rep(1:p, each = k)] * u[, rep(1:k, p)]), nrow = k) / m / EGu

  dlambda <- dlambda1 - n * dlambda2
  dlambda[lower.tri(dlambda)] <- 0
  dmu <- dmu1 + dmu2
  dAinv <- dAinv1 + dAinv2
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv))
}

#' @keywords internal
.fit.mstil.1 <- function(x, lambda, delta, Ainv, nu, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  if (!"maxitOptim" %in% names(control)) control$maxitOptim <- 10
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 0

  numGradSample <- control$numGradSample
  maxitOptim <- control$maxitOptim
  lambdaPenalty <- control$lambdaPenalty


  p <- ncol(lambda)
  k <- ncol(x)
  n <- nrow(x)
  u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)

  lik <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    return(sum(dmstil(x, lambda, delta, Ainv, nu, u, log.p = TRUE)) - lambdaPenalty * sum(abs(lambda)))
  }
  grad <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    grad <- .mstil.grad.1(x, lambda, delta, Ainv, nu, u)
    grad$dlambda <- grad$dlambda - lambdaPenalty * sign(lambda)
    gr <- c(as.vector(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)])
    return(gr)
  }
  param0 <- c(as.vector(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)])


  res <- stats::optim(param0, lik, grad, method = "BFGS", control = c(fnscale = -1, maxit = maxitOptim))
  param1 <- res$par
  lambda1 <- matrix(param1[1:(p * k)], nrow = k)
  delta1 <- param1[1:k + (p * k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (p * k) + 1):length(param1)]

  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu, value = res$value))
}

#' @keywords internal
.mstil.grad.2 <- function(x, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  if (!"finDiffStep" %in% names(control)) control$finDiffStep <- 1e-5
  numGradSample <- control$numGradSample
  finDiffStep <- control$finDiffStep
  n <- nrow(x)
  k <- ncol(x)
  if (missing(u)) u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  rsG <- rowSums(Gu)
  z <- t((t(x) - delta)) %*% t(Ainv)
  dnu2 <- (.dmvt2(u, rep(0, k), diag(k), nu = exp(log(nu) + finDiffStep)) - .dmvt2(u, rep(0, k), diag(k), nu = nu)) / finDiffStep
  dnu1 <- (.dmvt2(x, delta, Ainv, nu = exp(log(nu) + finDiffStep)) - .dmvt2(x, delta, Ainv, nu = nu)) / finDiffStep
  return(sum(dnu1) - n * sum(dnu2 * rsG) / sum(rsG))
}

#' @keywords internal
.fit.mstil.2 <- function(x, lambda, delta, Ainv, nu, control = list()) {
  if (!"stepSgd" %in% names(control)) control$stepSgd <- 1e-2
  if (!"iterSgd" %in% names(control)) control$iterSgd <- 1e2
  if (!"dimRateSgd" %in% names(control)) control$dimRateSgd <- 1e-2
  if (!"finDiffStep" %in% names(control)) control$finDiffStep <- 1e-5
  stepSgd <- control$stepSgd
  iterSgd <- control$iterSgd
  dimRateSgd <- control$dimRateSgd
  finDiffStep <- control$finDiffStep
  n <- nrow(x)
  k <- ncol(x)
  for (i in 1:iterSgd) {
    grad <- .mstil.grad.2(x, lambda, delta, Ainv, nu, control = control)
    if (is.na(grad)) grad <- 0
    nu <- exp(log(nu) + stepSgd * exp(-i * dimRateSgd) * grad / n / k)
    nu <- max(nu, finDiffStep)
  }
  return(nu)
}

#' @keywords internal
.mstil.grad.2.weighted <- function(x, w, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  if (!"finDiffStep" %in% names(control)) control$finDiffStep <- 1e-5
  numGradSample <- control$numGradSample
  finDiffStep <- control$finDiffStep
  k <- ncol(x)
  if (missing(u)) u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  rsG <- rowSums(Gu)
  dnu2 <- (.dmvt2(u, rep(0, k), diag(k), nu = exp(log(nu) + finDiffStep)) - .dmvt2(u, rep(0, k), diag(k), nu = nu)) / finDiffStep
  dnu1 <- (.dmvt2(x, delta, Ainv, nu = exp(log(nu) + finDiffStep)) - .dmvt2(x, delta, Ainv, nu = nu)) / finDiffStep
  return(sum(w * dnu1) - sum(w) * sum(dnu2 * rsG) / sum(rsG))
}

#' @keywords internal
.fit.mstil.2.weighted <- function(x, w, lambda, delta, Ainv, nu, control = list()) {
  if (!"stepSgd" %in% names(control)) control$stepSgd <- 1e-2
  if (!"iterSgd" %in% names(control)) control$iterSgd <- 1e2
  if (!"dimRateSgd" %in% names(control)) control$dimRateSgd <- 1e-2
  if (!"finDiffStep" %in% names(control)) control$finDiffStep <- 1e-5
  stepSgd <- control$stepSgd
  iterSgd <- control$iterSgd
  dimRateSgd <- control$dimRateSgd
  finDiffStep <- control$finDiffStep
  n <- nrow(x)
  k <- ncol(x)

  for (i in 1:iterSgd) {
    grad <- .mstil.grad.2.weighted(x, w, lambda, delta, Ainv, nu, control = control)
    if (is.na(grad)) grad <- 0
    nu <- exp(log(nu) + stepSgd * exp(-i * dimRateSgd) * grad / n / k)
    nu <- max(nu, finDiffStep)
  }
  return(nu)
}

#' @keywords internal
.mstil.grad.1.weighted <- function(x, w, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  numGradSample <- control$numGradSample


  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  if (missing(u)) u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)

  m <- nrow(u)
  y <- t(t(x) - delta)
  xA <- Ainv %*% t(y)
  z <- t(xA)
  xSx <- colSums(xA^2)
  C <- 1 / (1 + 1 / nu * xSx)

  lambdaz <- z %*% lambda
  lambdau <- u %*% lambda
  Gz <- stats::plogis(lambdaz, log.p = TRUE)
  Gu <- stats::plogis(lambdau, log.p = TRUE)
  gz <- stats::dlogis(lambdaz, log = TRUE)
  gu <- stats::dlogis(lambdau, log = TRUE)
  gGz <- exp(gz - Gz)
  gGu <- exp(gu - Gu)
  expGu <- exp(rowSums(Gu))
  EGu <- mean(expGu)

  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(w * (2 / nu * (y[, t2] * t(xA)[, t1]) * C))
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + sum(w) * diag(1 / diag(Ainv))
  dAinv1 <- grad
  dAinv2 <- t(t(matrix(colSums(w * y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)

  dmu1 <- colSums(w * t(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C))))
  dmu2 <- as.vector(-(t(colSums(w * gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda1 <- matrix(colSums(w * gGz[, rep(1:p, each = k)] * t(xA[rep(1:k, p), ])), nrow = k)
  dlambda2 <- matrix(colSums(expGu * gGu[, rep(1:p, each = k)] * u[, rep(1:k, p)]), nrow = k) / m / EGu

  dlambda <- dlambda1 - sum(w) * dlambda2
  dlambda[lower.tri(dlambda)] <- 0
  dmu <- dmu1 + dmu2
  dAinv <- dAinv1 + dAinv2
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv))
}

#' @keywords internal
.fit.mstil.1.weighted <- function(x, w, lambda, delta, Ainv, nu, control = list()) {
  if (!"numGradSample" %in% names(control)) control$numGradSample <- 1e4
  if (!"maxitOptim" %in% names(control)) control$maxitOptim <- 10
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 0


  numGradSample <- control$numGradSample
  maxitOptim <- control$maxitOptim
  lambdaPenalty <- control$lambdaPenalty
  p <- ncol(lambda)
  k <- ncol(x)
  n <- nrow(x)

  u <- mvtnorm::rmvt(numGradSample, delta = rep(0, k), sigma = diag(k), df = nu)
  lik <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    val <- sum(w * dmstil(x, lambda, delta, Ainv, nu, u, log.p = TRUE)) - lambdaPenalty * sum(abs(lambda))
    return(val)
  }
  grad <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    grad <- .mstil.grad.1.weighted(x, w, lambda, delta, Ainv, nu, u)
    grad$dlambda <- grad$dlambda - lambdaPenalty * sign(lambda)
    gr <- c(as.vector(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)])
    return(gr)
  }
  param0 <- c(as.vector(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)])

  res <- stats::optim(param0, lik, grad, method = "BFGS", control = c(fnscale = -1, maxit = maxitOptim))
  param1 <- res$par
  lambda1 <- matrix(param1[1:(p * k)], nrow = k)
  delta1 <- param1[1:k + (p * k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (p * k) + 1):length(param1)]

  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu, value = res$value))
}

#' @keywords internal
.mstil.r.grad <- function(x, lambda, delta, Ainv, nu) {
  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  y <- t(t(x) - delta)
  xA <- Ainv %*% t(y)
  z <- t(xA)
  xSx <- colSums(xA^2)
  C <- 1 / (1 + 1 / nu * xSx)

  lambdaz <- z %*% lambda
  Gz <- stats::plogis(lambdaz, log.p = TRUE)
  gz <- stats::dlogis(lambdaz, log = TRUE)

  gGz <- exp(gz - Gz)

  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(2 / nu * (y[, t2] * t(xA)[, t1]) * C)
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + n * diag(1 / diag(Ainv))
  dAinv1 <- grad
  dAinv2 <- t(t(matrix(colSums(y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)
  dmu1 <- rowSums(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C)))
  dmu2 <- as.vector(-(t(colSums(gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda1 <- diag(colSums(gGz * t(xA)))
  dnu1 <- sum((digamma((nu + k) / 2) - digamma(nu / 2)) / 2 - k / 2 / nu - 0.5 * (log1p(xSx / nu) + (nu + k) / (1 + xSx / nu) * xSx * (-1) / nu^2))
  dlambda <- dlambda1
  dmu <- dmu1 + dmu2
  dAinv <- dAinv1 + dAinv2
  dnu <- dnu1
  dlnu <- dnu1 * nu
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv, dnu = dnu, dlnu = dlnu))
}

#' @keywords internal
.mstil.r.grad.weighted <- function(x, w, lambda, delta, Ainv, nu) {
  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  y <- t(t(x) - delta)
  xA <- Ainv %*% t(y)
  z <- t(xA)
  xSx <- colSums(xA^2)
  C <- 1 / (1 + 1 / nu * xSx)
  lambdaz <- z %*% lambda
  Gz <- stats::plogis(lambdaz, log.p = TRUE)
  gz <- stats::dlogis(lambdaz, log = TRUE)
  gGz <- exp(gz - Gz)


  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(w * 2 / nu * (y[, t2] * t(xA)[, t1]) * C)
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + sum(w) * diag(1 / diag(Ainv))
  dAinv1 <- grad
  dAinv2 <- t(t(matrix(colSums(w * y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)

  dmu1 <- colSums(w * t(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C))))
  dmu2 <- as.vector(-(t(colSums(w * gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda1 <- diag(colSums(w * gGz * t(xA)))
  dnu1 <- sum(w * ((digamma((nu + k) / 2) - digamma(nu / 2)) / 2 - k / 2 / nu - 0.5 * (log1p(xSx / nu) + (nu + k) / (1 + xSx / nu) * xSx * (-1) / nu^2)))
  dlambda <- dlambda1
  dmu <- dmu1 + dmu2
  dAinv <- dAinv1 + dAinv2
  dnu <- dnu1
  dlnu <- dnu * nu
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv, dnu = dnu, dlnu = dlnu))
}

#' @keywords internal
.fit.fmmstil.r.weighted <- function(x, w, lambda, delta, Ainv, nu, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 0
  if (!"maxitOptimR" %in% names(control)) control$maxitOptimR <- 1e2

  maxitOptimR <- control$maxitOptimR
  lambdaPenalty <- control$lambdaPenalty

  k <- ncol(x)
  n <- nrow(x)
  lik <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    lnu <- param[length(param)]
    nu <- exp(lnu)
    res <- sum(w * dmstil.r(x, lambda, delta, Ainv, nu, log.p = TRUE) - lambdaPenalty * sum(abs(lambda)))
    return(res)
  }
  grad <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    lnu <- param[length(param)]
    nu <- exp(lnu)
    grad <- .mstil.r.grad.weighted(x, w, lambda, delta, Ainv, nu)
    grad$dlambda <- grad$dlambda - lambdaPenalty * sign(lambda)
    gr <- c(diag(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)], grad$dlnu)
    return(gr)
  }
  param0 <- c(diag(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)], log(nu))


  res <- stats::optim(param0, lik, grad, method = "BFGS", control = c(fnscale = -1, maxit = maxitOptimR))
  param1 <- res$par
  lambda1 <- diag(param1[1:k])
  delta1 <- param1[1:k + (k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (k) + 1):(length(param1) - 1)]
  lnu1 <- param1[length(param1)]
  nu1 <- exp(lnu1)
  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu1, value = res$value))
}

#' @keywords internal
.fmmstil.weight <- function(x, omega, lambda, delta, Ainv, nu, u, control = list()) {
  if (!"numLikSample" %in% names(control)) control$numLikSample <- 1e6
  numLikSample <- control$numLikSample
  K <- length(omega)
  k <- ncol(x)
  n <- nrow(x)
  w <- matrix(NA, nrow = n, ncol = K)
  
  if (missing(u)) {
    u <- list()
    for (i in 1:K) {
      u[[i]] <- mvtnorm::rmvt(numLikSample, delta = rep(0, k), sigma = diag(k), df = nu[[i]])
    }
  }
  for (i in 1:K) {
    w[, i] <- omega[[i]] * dmstil(x, lambda[[i]], delta[[i]], Ainv[[i]], nu[[i]], u[[i]], log.p = FALSE)
  }
  return(w)
}

#' @keywords internal
.fmmstil.r.weight <- function(x, omega, lambda, delta, Ainv, nu) {
  K <- length(omega)
  n <- nrow(x)
  w <- matrix(NA, nrow = n, ncol = K)

  for (i in 1:K) {
    w[, i] <- omega[[i]] * dmstil.r(x, lambda[[i]], delta[[i]], Ainv[[i]], nu[[i]], log.p = FALSE)
  }
  return(w)
}

#' @keywords internal
.default.init.param.method <- function(x) {
  param <- list(
    lambda = diag(stats::runif(ncol(x), -0.5, 0.5)),
    delta = x[sample(nrow(x), 1), ],
    Ainv = tryCatch(t(solve(chol(stats::cov(x)))),
      error = function(e) 1 / sqrt(stats::cov(x)),
      warning = function(w) 1 / sqrt(stats::cov(x))
    ),
    nu = 10
  )

  return(param)
}

#' @keywords internal
.default.init.cluster.method <- function(x, K) {
  prob <- stats::runif(K) + 1
  prob <- prob / sum(prob)
  init.cluster <- sample(K, nrow(x), replace = TRUE, prob = prob)
  return(init.cluster)
}

#' @keywords internal
.default.init.param.method2 <- function(x) {
  param <- list(
    lambda = diag(stats::runif(ncol(x), -0.5, 0.5)),
    delta = colMeans(x),
    Ainv = tryCatch(t(solve(chol(stats::cov(x)))),
                    error = function(e) 1 / sqrt(stats::cov(x)),
                    warning = function(w) 1 / sqrt(stats::cov(x))),
    nu = 10
  )
  return(param)
}

#' @keywords internal
.default.init.cluster.method2 <- function(x, K) {
  return(stats::kmeans(x, K)$cluster)
}

#' @keywords internal
.cluster.fmmstil.K.parallel.divisive <- function(x, K, ncore = 1, cluster0, init.cluster.method, init.param.method, show.progress = TRUE, control = list()) {
  if (!"lambdaPenalty" %in% names(control)) control$lambdaPenalty <- 1e-2
  if (!"cvgTolR" %in% names(control)) control$cvgTolR <- 1e-1
  if (ncore > parallel::detectCores()) {
    ncore <- parallel::detectCores()
    warning("Not enough available core")
  }
  if (missing(init.cluster.method)) init.cluster.method <- .default.init.cluster.method2
  if (missing(init.param.method)) init.param.method <- .default.init.param.method2
  
  if (show.progress) cat("\n", "MSTIL.R", "\t", "K : ", K, "\t", "Number of Trials : ", K - 1)
  seedFmmstil <- sample(.Machine$integer.max, 1)
  
  if (K == 2) cluster0 <- rep(1, nrow(x))
  ncore1 = min(ncore,(K-1))
  cluster0 <- factor(cluster0, labels = 1:K, levels = 1:K)
  
  initParamList <- list()
  for (trial in 1:(K - 1)){
    initParamList[[trial]] <- list(omega = list(), lambda = list(), delta = list(), Ainv = list(), nu = list())
    smallCluster <- factor(init.cluster.method(x[which(cluster0 == trial),], 2), labels = c(trial, K))
    initCluster <- cluster0
    initCluster[which(initCluster == trial)] <- smallCluster
    initCluster <- as.numeric(initCluster)
    
    initParamList[[trial]]$omega <- as.list(table(initCluster) / length(initCluster))
    for (kk in 1:K){
      initFit <- init.param.method(x[which(initCluster == unique(initCluster)[kk]),])
      initParamList[[trial]]$lambda[[kk]] = initFit$lambda
      initParamList[[trial]]$delta[[kk]] = initFit$delta
      initParamList[[trial]]$Ainv[[kk]] = initFit$Ainv
      initParamList[[trial]]$nu[[kk]] = initFit$nu
    }
  }
  
  
  fit.fmmstil.r.seed.windows <- function(trial, data, K, initParamList, control = list()){
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
  resRec <- parallel::parLapply(NULL, 1:(K-1), fit.fmmstil.r.seed.windows, control = control, data = x, K = K, initParamList = initParamList)
  parallel::stopCluster(cl)
  set.seed(seedFmmstil)
  
  
  ICLRec <- c()
  for (i in 1:(K-1)) ICLRec <- c(ICLRec, resRec[[i]]$ICL)
  res1Best <- resRec[[which.max(ICLRec)]]
  par <- res1Best$par[[which.max(res1Best$logLik)]]
  
  if (show.progress) cat("\t", "Max ICL : ", (round(max(ICLRec), 2)))
  
  
  if (show.progress) cat("\n", "MSTIL  ", "\t", "K : ", K, "\t")
  res2Best <- tryCatch(fit.fmmstil.parallel(x, K, ncore = ncore, param = par, show.progress = FALSE, control = control),
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



