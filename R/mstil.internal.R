#' @keywords internal
dmvt2 <- function(x, delta, Ainv, df = 1, log = TRUE) {
  p <- ncol(x)
  R.x_m <- Ainv %*% (t(x) - delta)
  rss <- colSums(R.x_m^2)
  logretval <- lgamma((p + df) / 2) - (lgamma(df / 2) - sum(log(abs(diag(Ainv)))) +
    p / 2 * log(pi * df)) - 0.5 * (df + p) * log1p(rss / df)
  if (log) {
    logretval
  } else {
    exp(logretval)
  }
}
#' @keywords internal
dmstil.grad <- function(x, lambda, delta, Ainv, nu, u, sample.mc = 10000) {
  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }
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
  exp_Gu <- exp(rowSums(Gu))
  E_Gu <- mean(exp_Gu)

  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(2 / nu * (y[, t2] * t(xA)[, t1]) * C)
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + n * diag(1 / diag(Ainv))
  dAinv_1 <- grad
  dAinv_2 <- t(t(matrix(colSums(y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)
  dmu_1 <- rowSums(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C)))
  dmu_2 <- as.vector(-(t(colSums(gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda_1 <- matrix(colSums(gGz[, rep(1:p, each = k)] * t(xA[rep(1:k, p), ])), nrow = k)
  dlambda_2 <- matrix(colSums(exp_Gu * gGu[, rep(1:p, each = k)] * u[, rep(1:k, p)]), nrow = k) / m / E_Gu

  dlambda <- dlambda_1 - n * dlambda_2
  dlambda[lower.tri(dlambda)] <- 0
  dmu <- dmu_1 + dmu_2
  dAinv <- dAinv_1 + dAinv_2
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv))
}
#' @keywords internal
fit.mstil.qcmle <- function(x, lambda, delta, Ainv, nu, u, sample.mc = 10000, maxit = 10, lambda.penalty = 0, max.lambda = Inf) {
  p <- ncol(lambda)
  k <- ncol(x)
  n <- nrow(x)
  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }

  lik <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    return(sum(dmstil(x, lambda, delta, Ainv, nu, u, log.p = TRUE)) - sum(exp(lambda.penalty * lambda^2)))
  }
  grad <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    grad <- dmstil.grad(x, lambda, delta, Ainv, nu, u)
    grad$lambda <- grad$lambda - lambda.penalty * 2 * lambda * exp(lambda.penalty * lambda^2)
    gr <- c(as.vector(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)])
    return(gr)
  }
  param0 <- c(as.vector(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)])

  if (is.infinite(max.lambda)) optim.method = 'BFGS'
  else optim.method = 'L-BFGS-B'
  lower = c(rep(-max.lambda, p * k), rep(-Inf, length(param0) - p * k))
  upper = -lower
  res <- stats::optim(param0, lik, grad, method = optim.method, control = c(fnscale = -1, maxit = maxit), lower = lower, upper = upper )
  param1 <- res$par
  lambda1 <- matrix(param1[1:(p * k)], nrow = k)
  delta1 <- param1[1:k + (p * k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (p * k) + 1):length(param1)]

  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu, value = res$value))
}
#' @keywords internal
dmstil.grad.nu <- function(x, lambda, delta, Ainv, nu, u, sample.mc = 10000, tol = 1e-5) {
  n <- nrow(x)
  k <- ncol(x)
  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  rsG <- rowSums(Gu)
  z <- t((t(x) - delta)) %*% t(Ainv)
  dnu_u <- (dmvt2(u, rep(0, k), diag(k), df = exp(log(nu) + tol)) - dmvt2(u, rep(0, k), diag(k), df = nu)) / tol
  dnu_x <- (dmvt2(x, delta, Ainv, df = exp(log(nu) + tol)) - dmvt2(x, delta, Ainv, df = nu)) / tol
  return(sum(dnu_x) - n * sum(dnu_u * rsG) / sum(rsG))
}
#' @keywords internal
fit.mstil.sgd <- function(x, lambda, delta, Ainv, nu, tol = 1e-5, step.size = 0.1, dim.rate = 0.1, iter.sgd = 100, sample.mc = 10000) {
  for (i in 1:iter.sgd) {
    grad <- dmstil.grad.nu(x, lambda, delta, Ainv, nu, sample.mc = sample.mc)
    if (is.na(grad)) grad <- 0
    nu <- exp(log(nu) + step.size * exp(-i * dim.rate) * grad)
    nu <- max(nu, tol)
  }
  return(nu)
}
#' @keywords internal
dmstil.grad.weighted.nu <- function(x, w, lambda, delta, Ainv, nu, u, sample.mc = 10000, tol = 1e-5) {
  k <- ncol(x)
  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }
  Gu <- stats::plogis(u %*% lambda, log.p = TRUE)
  rsG <- rowSums(Gu)
  dnu_u <- (dmvt2(u, rep(0, k), diag(k), df = exp(log(nu) + tol)) - dmvt2(u, rep(0, k), diag(k), df = nu)) / tol
  dnu_x <- (dmvt2(x, delta, Ainv, df = exp(log(nu) + tol)) - dmvt2(x, delta, Ainv, df = nu)) / tol
  return(sum(w * dnu_x) - sum(w) * sum(dnu_u * rsG) / sum(rsG))
}
#' @keywords internal
fit.mstil.weighted.sgd <- function(x, w, lambda, delta, Ainv, nu, tol = 1e-5, step.size = 0.1, dim.rate = 0.1, iter.sgd = 100, sample.mc = 10000) {
  for (i in 1:iter.sgd) {
    grad <- dmstil.grad.weighted.nu(x, w, lambda, delta, Ainv, nu, sample.mc = sample.mc)
    if (is.na(grad)) grad <- 0
    nu <- exp(log(nu) + step.size * exp(-i * dim.rate) * grad)
    nu <- max(nu, tol)
  }
  return(nu)
}
#' @keywords internal
dmstil.grad.weighted <- function(x, w, lambda, delta, Ainv, nu, u, sample.mc = 10000) {
  n <- nrow(x)
  k <- ncol(x)
  p <- ncol(lambda)
  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }
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
  exp_Gu <- exp(rowSums(Gu))
  E_Gu <- mean(exp_Gu)

  t1 <- matrix(rep(1:k, k), nrow = k)[lower.tri(diag(k), diag = TRUE)]
  t2 <- rep(1:k, k:1)
  grad <- matrix(0, nrow = k, ncol = k)
  grad_ <- colSums(w * (2 / nu * (y[, t2] * t(xA)[, t1]) * C))
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + sum(w) * diag(1 / diag(Ainv))
  dAinv_1 <- grad
  dAinv_2 <- t(t(matrix(colSums(w * y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)

  dmu_1 <- colSums(w * t(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C))))
  dmu_2 <- as.vector(-(t(colSums(w * gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda_1 <- matrix(colSums(w * gGz[, rep(1:p, each = k)] * t(xA[rep(1:k, p), ])), nrow = k)
  dlambda_2 <- matrix(colSums(exp_Gu * gGu[, rep(1:p, each = k)] * u[, rep(1:k, p)]), nrow = k) / m / E_Gu

  dlambda <- dlambda_1 - sum(w) * dlambda_2
  dlambda[lower.tri(dlambda)] <- 0
  dmu <- dmu_1 + dmu_2
  dAinv <- dAinv_1 + dAinv_2
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv))
}
#' @keywords internal
fit.mstil.weighted.qmle <- function(x, w, lambda, delta, Ainv, nu, u, sample.mc = 10000, maxit = 10, lambda.penalty = 0, max.lambda = Inf) {
  p <- ncol(lambda)
  k <- ncol(x)
  n <- nrow(x)

  if (missing(u)) {
    u <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu)
  }
  lik <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    return(sum(w * dmstil(x, lambda, delta, Ainv, nu, u, log.p = TRUE)) - sum(exp(lambda.penalty * lambda^2)))
  }
  grad <- function(param) {
    lambda <- matrix(param[1:(p * k)], nrow = k)
    delta <- param[1:k + (p * k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + (p * k) + 1):length(param)]
    grad <- dmstil.grad.weighted(x, w, lambda, delta, Ainv, nu, u)
    grad$lambda <- grad$lambda - lambda.penalty * 2 * lambda * exp(lambda.penalty * lambda^2)
    gr <- c(as.vector(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)])
    return(gr)
  }
  param0 <- c(as.vector(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)])
  
  if (is.infinite(max.lambda)) optim.method = 'BFGS'
  else optim.method = 'L-BFGS-B'
  lower = c(rep(-max.lambda, p * k), rep(-Inf, length(param0) - p * k))
  upper = -lower
  res <- stats::optim(param0, lik, grad, method = optim.method, control = c(fnscale = -1, maxit = maxit), lower = lower, upper = upper)
  param1 <- res$par
  lambda1 <- matrix(param1[1:(p * k)], nrow = k)
  delta1 <- param1[1:k + (p * k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (p * k) + 1):length(param1)]

  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu, value = res$value))
}
#' @keywords internal
mstil.weight <- function(x, omega, lambda, delta, Ainv, nu, u, sample.mc = 1000000) {
  K <- length(omega)
  k <- ncol(x)
  n <- nrow(x)
  w <- matrix(NA, nrow = n, ncol = K)

  if (missing(u)) {
    u <- list()
    for (i in 1:K) u[[i]] <- mvtnorm::rmvt(sample.mc, delta = rep(0, k), sigma = diag(k), df = nu[[i]])
  }
  for (i in 1:K) {
    w[, i] <- omega[i] * dmstil(x, lambda[[i]], delta[[i]], Ainv[[i]], nu[[i]], u[[i]], log.p = FALSE)
  }
  return(w)
}
#' @keywords internal
dmstil.r.grad <- function(x, lambda, delta, Ainv, nu) {
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
  dAinv_1 <- grad
  dAinv_2 <- t(t(matrix(colSums(y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)
  dmu_1 <- rowSums(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C)))
  dmu_2 <- as.vector(-(t(colSums(gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda_1 <- diag(colSums(gGz * t(xA)))
  dnu_1 <- sum((digamma((nu + k) / 2) - digamma(nu / 2)) / 2 - k / 2 / nu - 0.5 * (log1p(xSx / nu) + (nu + k) / (1 + xSx / nu) * xSx * (-1) / nu^2))
  dlambda <- dlambda_1
  dmu <- dmu_1 + dmu_2
  dAinv <- dAinv_1 + dAinv_2
  dnu <- dnu_1
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv, dnu = dnu))
}
#' @keywords internal
fit.mstil.r.weighted <- function(x, w, lambda, delta, Ainv, nu, maxit = 1000, max.lambda = Inf) {
  k <- ncol(x)
  n <- nrow(x)
  lik <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    nu <- param[length(param)]
    nu <- max(1, nu)
    res <- sum(w * dmstil.r(x, lambda, delta, Ainv, nu, log.p = TRUE))
    return(res)
  }
  grad <- function(param) {
    lambda <- diag(param[1:k])
    delta <- param[1:k + (k)]
    Ainv <- matrix(0, nrow = k, ncol = k)
    Ainv[lower.tri(Ainv, diag = TRUE)] <- param[(k + k + 1):(length(param) - 1)]
    nu <- param[length(param)]
    grad <- dmstil.r.grad.weighted(x, w, lambda, delta, Ainv, nu)
    gr <- c(diag(grad$dlambda), grad$dmu, grad$dAinv[lower.tri(diag(k), diag = TRUE)], grad$dnu)
    return(gr)
  }
  param0 <- c(diag(lambda), delta, Ainv[lower.tri(diag(k), diag = TRUE)], nu)
  
  if (is.infinite(max.lambda)) optim.method = 'BFGS'
  else optim.method = 'L-BFGS-B'
  lower = c(rep(-max.lambda, k), rep(-Inf, length(param0) - k))
  upper = -lower
  
  res <- stats::optim(param0, lik, grad, method = optim.method, control = c(fnscale = -1, maxit = maxit), lower = lower, upper = upper)
  param1 <- res$par
  lambda1 <- diag(param1[1:k])
  delta1 <- param1[1:k + (k)]
  Ainv1 <- matrix(0, nrow = k, ncol = k)
  Ainv1[lower.tri(Ainv1, diag = TRUE)] <- param1[(k + (k) + 1):(length(param1) - 1)]
  nu1 <- param1[length(param1)]

  return(list(lambda = lambda1, delta = delta1, Ainv = Ainv1, nu = nu1, value = res$value))
}
#' @keywords internal
dmstil.r.grad.weighted <- function(x, w, lambda, delta, Ainv, nu) {
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
  grad_ <- colSums( w * 2 / nu * (y[, t2] * t(xA)[, t1]) * C)
  grad[lower.tri(grad, diag = TRUE)] <- grad_
  grad <- -(nu + k) / 2 * grad + sum(w) * diag(1 / diag(Ainv))
  dAinv_1 <- grad
  dAinv_2 <- t(t(matrix(colSums(w * y[, rep(1:k, each = p)] * gGz[, rep(1:p, k)]), nrow = p)) %*% t(lambda)) * lower.tri(diag(k), diag = TRUE)
  
  dmu_1 <- colSums(w * t(-(nu + k) / 2 / nu * (-2) * (crossprod(Ainv) %*% t(y * C))))
  dmu_2 <- as.vector(-(t(colSums(w * gGz)) %*% (t(lambda) %*% Ainv)))
  dlambda_1 <- diag(colSums(w * gGz * t(xA)))
  dnu_1 <- sum(w * ((digamma((nu + k) / 2) - digamma(nu / 2)) / 2 - k / 2 / nu - 0.5 * (log1p(xSx / nu) + (nu + k) / (1 + xSx / nu) * xSx * (-1) / nu^2)))
  dlambda <- dlambda_1
  dmu <- dmu_1 + dmu_2
  dAinv <- dAinv_1 + dAinv_2
  dnu <- dnu_1
  return(list(dlambda = dlambda, dmu = dmu, dAinv = dAinv, dnu = dnu))
}
#' @keywords internal
mstil.r.weight <- function(x, omega, lambda, delta, Ainv, nu) {
  K <- length(omega)
  n <- nrow(x)
  w <- matrix(NA, nrow = n, ncol = K)

  for (i in 1:K) {
    w[, i] <- omega[i] * dmstil.r(x, lambda[[i]], delta[[i]], Ainv[[i]], nu[[i]], log.p = FALSE)
  }
  return(w)
}
