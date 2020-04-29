generate_skew = function( p, data_size ){
  m1 = matrix(runif(p*p,-1,1),nrow=p)*3
  
  S = LaplacesDemon::rinvwishart(2*p+1,diag(2*p))
  a = mvtnorm::rmvt( data_size, delta = runif(2*p, -3,3), sigma = S, df = runif(1, 5, 20) )
  d = a[,1:p]+abs(a[,(p+1):(2*p)])%*%m1
  return(d)
}
set.seed(3)
K = 2
k = 5
membership = c()
n1 = round(runif(1, 100,500))*5
membership = c(membership,rep(1,n1))
dat = generate_skew(k, n1)
for ( i in 2:K){
  n1 = round(runif(1, 100,500))*5
  membership = c(membership,rep(i,n1))
  dat_ = generate_skew(k, n1)
  dat = rbind(dat, dat_)
}
set.seed(NULL)

plot(0,cex=0.1,xlim=range(dat[,1]),ylim=range(dat[,2]))
for (i in 1: K){
  lines(dat[which(membership==i),],col=(i+1),type='p',cex=0.5)
}

res1 = mixtools::mvnormalmixEM(dat, k = K)
print('logL of GMM')
print(res1$loglik)


lambda = list()
delta = list()
Ainv = list()
nu = list()
for ( i in 1:K){
    res_ = fit.mstil(dat[which(membership==i),],print.progress = FALSE)
  mm = which.max(res_$logL)
  lambda[[i]] = res_$par[[mm]]$lambda
  delta[[i]] = res_$par[[mm]]$delta
  Ainv[[i]] = res_$par[[mm]]$Ainv
  nu[[i]] = res_$par[[mm]]$nu
}
omega = table(membership)/length(membership)

print('given true membership')
res_true = fit.fmmstil(dat, K, omega, lambda, delta, Ainv, nu)

print('restricted fmmstil')
res.fmmstil.r = fit.fmmstil.r(dat,K)

print('fmmstil')
res.fmmstil = fit.fmmstil(dat,K,
                          res.fmmstil.r$par[[length(res.fmmstil.r$par)]]$omega,
                          res.fmmstil.r$par[[length(res.fmmstil.r$par)]]$lambda,
                          res.fmmstil.r$par[[length(res.fmmstil.r$par)]]$delta,
                          res.fmmstil.r$par[[length(res.fmmstil.r$par)]]$Ainv,
                          res.fmmstil.r$par[[length(res.fmmstil.r$par)]]$nu,
                          maxit.bfgs = 10,maxit.qmle = 1,convergence.n = 5)

print('fmcfust')
res.cfust = EMMIXcskew::fmcfust(K,dat)
print('fmmst')
res.mmst = EMMIXcskew::fmmst(K,dat)