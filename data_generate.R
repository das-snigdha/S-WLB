# Code to generate data under a survey design with unequal probability of sampling

# required package
library(mvtnorm)

# Simulation 1: Normal model
# Generate the variable of interest X by using selection variable Z, where
# (X, Z) jointly follow a Bivariate normal distribution
# Inclusion probabilities are calculated using a probit link function
# mu_x, s_x : mean and standard deviation of X 
# mu_z, s_z : mean and standard deviation of Z
# rho : correlation between X and Z
# b0, b1: parameters used to calculate inclusion probabilities
# n : sample size
# N : population size
data_normal = function(mu_x, mu_z, s_x, s_z, rho, b0, b1, n, N){
  
  # mean vector
  mu = c(mu_x, mu_z)
  
  # covariance matrix
  S = matrix(c(s_x^2, rho*s_x*s_z, rho*s_x*s_z, s_z^2), nrow = 2)
  
  # generate the population
  data.pop = rmvnorm(n = N, mean = mu, sigma = S)
  colnames(data.pop) = c("X", "Z")
  
  # selection probabilities using probit link
  # pi = \Phi(b_0 + b_1*Z)
  sel.prob = pnorm(b0 + b1*data.pop[ , "Z"])
  
  # sample indicators
  sel.ind = sample(1:N, size = n, prob = sel.prob)
  wts = 1/sel.prob
  
  # weights corresponding to the selected sample (adding up to N)
  sample.wts = wts[sel.ind]
  sample.wts = N * sample.wts/sum(sample.wts)
  
  # values of the selected sample
  x = data.pop[sel.ind , "X"]
  
  df.pop = data.frame("X" = data.pop[ , "X"], "Z" = data.pop[ , "Z"], 
                      "sel.prob" = sel.prob)
  
  out = list("x" = x, "w" = sample.wts, "population" = df.pop)
  return(out)
}


# Simulation 2: Probit Regression Model
# Generate the bivariate data (X,Y) and selection variable Z, as follows:
# Covariate X follows a Normal distribution
# (V, Z) jointly follow a Bivariate normal distribution with mean of V being X\beta
# Response Y = 1 if V > 0, 0 if V < 0
# Inclusion probabilities are calculated using a probit link function
# beta_reg : regression parameter
# mu_x, s_x : mean and standard deviation of X 
# mu_z, s_z : mean and standard deviation of Z
# s_v : standard deviation of V
# rho : correlation between V and Z
# b0, b1: parameters used to calculate inclusion probabilities
# n : sample size
# N : population size
data_probit = function(beta_reg, mu_x = 1, s_x = 0.1, mu_z = 0, s_z = 1, s_v = 1, 
                       rho, b0, b1, n, N){
  
  X = rnorm(N, mean = mu_x, sd = s_x)
  mu_v = X*beta_reg
  
  # covariance matrix
  S = matrix(c(s_v^2, rho*s_v*s_z, rho*s_v*s_z, s_z^2), nrow = 2)
  
  # generate the population
  data.pop = rmvnorm(n = N, mean = c(0, 0), sigma = S)
  colnames(data.pop) = c("V", "Z")
  data.pop[, "V"] = mu_v + data.pop[, "V"] 
  data.pop[, "Z"] = mu_z + data.pop[, "Z"] 
  
  Y = ifelse(data.pop[ , "V"] >= 0, 1, 0)
  
  # selection probabilities using probit link
  # pi = \Phi(b_0 + b_1*Z)
  sel.prob = pnorm(b0 + b1*data.pop[ , "Z"])
  
  # sample indicators
  sel.ind = sample(1:N, size = n, prob = sel.prob)
  wts = 1/sel.prob
  
  # weights corresponding to the selected sample (adding up to N)
  sample.wts = wts[sel.ind]
  sample.wts = N * sample.wts/sum(sample.wts)
  
  # response of the selected sample
  y_s = Y[sel.ind]
  x_s = X[sel.ind]
  
  df.pop = data.frame("Y" = Y, "X" = X, "V" = data.pop[ , "V"],
                      "Z" = data.pop[ , "Z"], "sel.prob" = sel.prob)
  
  out = list("y" = y_s, "x" = x_s, "w" = sample.wts, "population" = df.pop)
  return(out)
}