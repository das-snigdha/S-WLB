# Estimators for Simulation 1: Normal model

# required package
library(extraDistr)

# Weighted Finite Population Bayesian Bootstrap (Cohen, 1997)
# y : observed sample of size n
# w : corresponding sampling weights
# N : population size
WFPBB.normal = function(y, w, N){
  
  n = length(y) ; m = N - n
  l = rep(0, n); y.star = rep(NA, m)
  
  for(k in 1:(N-n)){
    denom = m*(1 + (k-1)/n)
    
    w.star = ((w - 1) + l * (m/n))/ denom
    
    y.star[k] = sample(y, size = 1, prob = w.star)
    
    ind = which(y.star[k] == y)
    l[ind] = l[ind] + 1
  }
  
  y.all = c(y, y.star)
  
  y.rep = sample(y.all, size = n)
  
  return(y.rep)
}

# Pseudo Maximum Likelihood Estimator
# y : observed sample of size n
# w : corresponding sampling weights
PMLE.normal = function(y, w){
  
  n = length(y)
  w.tilde = n * (w/sum(w))
  
  mu = sum(w.tilde * y)/n
  
  sigma.sq = sum( (w.tilde^2) * ((y - mu)^2) ) / (n^2)
  
  CI = c(-1.96, 1.96)*sqrt(sigma.sq)
  CI = mu + CI
  
  return(list(mu.PMLE = mu, sigma.sq.PMLE = sigma.sq, CI.PMLE = CI))
}

# Unweighted Bayesian Estimator
# y : observed sample of size n
UBE.normal = function(y){
  
  n = length(y)
  mu = mean(y)
  s.sq = sum((y - mu)^2)/(n-1)
  
  sigma.sq = ((n-1)*s.sq)/(n*(n-3))
  
  CI = c(-1.96, 1.96)*sqrt(sigma.sq)
  CI = mu + CI
  
  return(list(mu.UBE = mu, sigma.sq.UBE = sigma.sq, CI.UBE = CI))
}

# Bayesian Pseudo Posterior Estimator
# y : observed sample of size n
# w : corresponding sampling weights
BPPE.normal = function(y, w){
  
  n = length(y)
  w.tilde = n * (w/sum(w))
  
  mu = sum(w.tilde * y)/n
  s.sq = sum( w.tilde * ((y - mu)^2) ) /(n-1)
  
  sigma.sq = ((n-1)*s.sq) / (n*(n-3))
  
  CI = c(-1.96, 1.96)*sqrt(sigma.sq)
  CI = mu + CI
  
  return(list(mu.BPPE = mu, sigma.sq.BPPE = sigma.sq, CI.BPPE = CI))
  
}

# Estimator from Weighted Bayesian Bootstrap
# y : observed sample of size n
# w : corresponding sampling weights
# N : population size
# M : number of posterior draws
WBB.normal = function(y, w, N, M){
  
  n = length(y)
  zbar = rep(NA, M)
  sigma.sq.y = rep(NA, M)
  
  # Algorithm 2 from Gunawan et al. (2020)
  for(i in 1:M){
    z = WFPBB.normal(y = y, w = w, N = N)
    
    zbar[i] = mean(z)
    sz.sq = sum((z - zbar[i])^2)
    
    sigma.sq.y[i] = 1/rgamma(1, shape = 0.5*(n-1), rate = sz.sq*0.5)
    
    if(i %% 50 == 0){
      print(paste("iteration :", i))
    }
    
  }
  mu = mean(zbar)
  sigma.sq = mean(sigma.sq.y)/n + mean((zbar - mu)^2)
  
  CI = c(-1.96, 1.96)*sqrt(sigma.sq)
  CI = mu + CI
  
  return(list(mu.WBB = mu, sigma.sq.WBB = sigma.sq, CI.WBB = CI))
  
}

# Survey adjusted Weighted Likelihood Bootstrap
# y : observed sample of size n
# w : corresponding sampling weights
# B : number of Bootstrap samples
SWLB.normal = function(y, w, B){
  n = length(y)
  w.tilde = n* (w/sum(w))
  
  mu.samp = rep(NA, B)
  
  for(i in 1:B){
    
    g = rgamma(n, shape = 1, rate = (1/w.tilde))
    g = g/sum(g)
    
    mu.samp[i] = sum(g*y)
  }
  
  mu = mean(mu.samp)
  sigma.sq = var(mu.samp)
  
  CI = quantile(mu.samp, probs = c(0.025,0.975))
  
  return(list(mu.SWLB = mu, sigma.sq.SWLB = sigma.sq, CI.SWLB = CI))
}