# Estimators for Simulation 2: Probit regression model

# required packages
library(extraDistr)
library(rstanarm)
library(sandwich)

# Weighted Finite Population Bayesian Bootstrap (Cohen, 1997)
# x, y : observed sample of size n
# w : corresponding sampling weights
# N : population size
WFPBB.probit = function(y, x, w, N){
  
  n = length(y) ; m = N - n
  l = rep(0, n); y.star = x.star = rep(NA, m)
  
  for(k in 1:(N-n)){
    denom = m*(1 + (k-1)/n)
    
    w.star = ((w - 1) + l * (m/n))/ denom
    
    ind = sample(1:n, size = 1, prob = w.star)
    
    y.star[k] = y[ind]
    x.star[k] = x[ind]
    l[ind] = l[ind] + 1
  }
  
  y.all = c(y, y.star)
  x.all = c(x, x.star)
  
  samp.ind = sample(1:N, size = n)
  
  y.rep = y.all[samp.ind]
  x.rep = x.all[samp.ind]
  
  df = data.frame(x.rep = x.rep, y.rep = y.rep)
  
  return(df)
}

# Pseudo Maximum Likelihood Estimator
# x, y : observed sample of size n
# w : corresponding sampling weights
PMLE.probit = function(y, x, w){
  
  n = length(y)
  w.tilde = n * (w/sum(w))
  
  suppressMessages(suppressWarnings({
    capture.output({
      fit = glm(y ~ x-1, family = quasibinomial(link = "probit"),
                weights = w.tilde, maxit = 200, epsilon = 1e-5)
    })
  }))
  
  
  beta.fit = fit$coefficients
  sigma.sq = c(sandwich(fit))
  
  CI = c(-1.96, 1.96)*sqrt(sigma.sq)
  CI = beta.fit + CI
  
  return(list(beta.PMLE = beta.fit, sigma.sq.PMLE = sigma.sq, CI.PMLE = CI))
}

# Unweighted Bayesian Estimator
# x, y : observed sample of size n
# B : number of MCMC samples required
# Burn.in : Burn in period of the MCMC chain
UBE.probit = function(y, x, B, Burn.in){
  
  df = data.frame(y = y, x = x)
  
  
  suppressMessages(suppressWarnings({
    capture.output({
      fit = stan_glm(y ~ x-1, data = df, family = binomial(link = "probit"),
                     init = 0, algorithm = "sampling",
                     prior = default_prior_coef(family), chains = 1, 
                     warmup = Burn.in, iter = B+Burn.in)
    })
  }))
  
  beta.fit = c(as.matrix(fit))
  beta.hat = mean(beta.fit)
  
  sigma.sq = var(beta.fit)
  
  CI = quantile(beta.fit, probs = c(0.025, 0.975))
  
  return(list(beta.UBE = beta.hat, sigma.sq.UBE = sigma.sq, CI.UBE = CI))
}

# Bayesian Pseudo Posterior Estimator
# x, y : observed sample of size n
# w : corresponding sampling weights
# B : number of MCMC samples required
# Burn.in : Burn in period of the MCMC chain
BPPE.probit = function(y, x, w, B, Burn.in){
  
  n = length(y)
  w.tilde = n * (w/sum(w))
  df = data.frame(y = y, x = x)
  
  suppressMessages(suppressWarnings({
    capture.output({
      fit = stan_glm(y ~ x-1, data = df, family = binomial(link = "probit"),
                     weights = w.tilde, init = 0, algorithm = "sampling",
                     prior = default_prior_coef(family), chains = 1, 
                     warmup = Burn.in, iter = B+Burn.in)
    })
  }))
  
  beta.fit = c(as.matrix(fit))
  beta.hat = mean(beta.fit)
  
  sigma.sq = var(beta.fit)
  
  CI = quantile(beta.fit, probs = c(0.025, 0.975))
  
  return(list(beta.BPPE = beta.hat, sigma.sq.BPPE = sigma.sq, CI.BPPE = CI))
  
}

# Estimator from Weighted Bayesian Bootstrap
# x, y : observed sample of size n
# w : corresponding sampling weights
# N : population size
# J : number of Bootstrap replicates
# M : number of MCMC samples required
# Burn.in : Burn in period of the MCMC chain
WBB.probit = function(y, x, w, N, J, M, Burn.in){
  
  n = length(y)
  beta.mat = matrix(NA, nrow = J, ncol = M)
  beta.fit = rep(NA, J*M)
  
  # Algorithm 3 from Gunawan et al. (2020)
  for(j in 1:J){
    z = WFPBB.probit(y = y, x = x, w = w, N = N)
    
    suppressMessages(suppressWarnings({
      capture.output({
        fit = stan_glm(y.rep ~ x.rep-1, data = z, family = binomial(link = "probit"),
                       init = 0, algorithm = "sampling",
                       prior = default_prior_coef(family), chains = 1, 
                       warmup = Burn.in, iter = M+Burn.in)
      })
    }))
    
    beta.mat[j, ] = c(as.matrix(fit))
    
    if(j %% 10 == 0){
      print(paste("iteration :", j))
    }
  }
  
  beta.fit = c(beta.mat)
  beta.hat = mean(beta.fit)
  
  sigma.sq = var(beta.fit)
  
  CI = quantile(beta.fit, probs = c(0.025, 0.975))
  
  return(list(beta.WBB = beta.hat, sigma.sq.WBB = sigma.sq, CI.WBB = CI))
  
}

# Survey adjusted Weighted Likelihood Bootstrap
# x, y : observed sample of size n
# w : corresponding sampling weights
# B : number of Bootstrap samples
SWLB.probit = function(y, x, w, B){
  n = length(y)
  w.tilde = n* (w/sum(w))
  
  beta.fit = rep(NA, B)
  
  for(i in 1:B){
    
    g = rgamma(n, shape = 1, rate = (1/w.tilde))
    g = g/sum(g)
    
    df = data.frame(y = y, x = x, w = g)
    
    fit = glm(y ~ x-1, weights = g, family = quasibinomial(link = "probit"),
              maxit = 200, epsilon = 1e-5)

    beta.fit[i] = fit$coefficients
  }
  
  beta.hat = mean(beta.fit)
  sigma.sq = var(beta.fit)
  
  CI = quantile(beta.fit, probs = c(0.025, 0.975))
  
  return(list(beta.SWLB = beta.hat, sigma.sq.SWLB = sigma.sq, CI.SWLB = CI))
}