Scalable Efficient Inference in Complex Surveys through Targeted
Resampling of Weights
================

Implementation of modeling framework proposed in the paper [Das, S.,
Bandyopadhyay, D., and Pati, D. (2025+) Scalable Efficient Inference in
Complex Surveys through Targeted
Resampling](https://arxiv.org/abs/2504.11636).

## Overview

Survey data often arises from complex sampling designs with unequal
inclusion probabilities. When sampling is informative, traditional
inference methods yield biased estimators and poor coverage. For
efficient inference on survey data, we propose the *Survey-adjusted
Weighted Likelihood Bootstrap* (S-WLB), which resamples weights from a
from a carefully chosen distribution centered around the underlying
sampling weights. S-WLB is computationally efficient, theoretically
consistent, and delivers finite-sample uncertainty intervals which are
proven to be asymptotically valid.

## Functions

- **normal_model.R** - Functions to get the estimators when data is
  generated from a normal model. Estimates are obtained from the
  following functions:

  - `SWLB.normal` - Our proposed Survey-adjusted Weighted Likelihood
    Bootstrap.
  - `PMLE.normal` - Pseudo Maximum Likelihood Estimator.
  - `WBB.normal` - Weighted Bayesian Bootstrap, Algorithm 2 of [Gunawan
    et
    al. (2020)](https://onlinelibrary.wiley.com/doi/abs/10.1111/anzs.12284).
  - `BPPE.normal` - Bayesian Pseudo Posterior Estimator ([Savitsky &
    Toth,
    2016](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Bayesian-estimation-under-informative-sampling/10.1214/16-EJS1153.full)).
  - `UBE.normal` - Unweighted Bayesian Estimator.

- **probit_model.R** - Functions to get the estimators when data is
  generated from a probit regression model. Estimates are obtained from
  the following functions:

  - `SWLB.probit` - Our proposed Survey-adjusted Weighted Likelihood
    Bootstrap.
  - `PMLE.probit` - Pseudo Maximum Likelihood Estimator.
  - `WBB.probit` - Weighted Bayesian Bootstrap, Algorithm 3 of [Gunawan
    et
    al. (2020)](https://onlinelibrary.wiley.com/doi/abs/10.1111/anzs.12284).
  - `BPPE.probit` - Bayesian Pseudo Posterior Estimator ([Savitsky &
    Toth,
    2016](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Bayesian-estimation-under-informative-sampling/10.1214/16-EJS1153.full)).
  - `UBE.probit` - Unweighted Bayesian Estimator.

- **data_generate.R** - Functions to generate data under a survey design
  with unequal probability of sampling. The functions `data_normal` and
  `data_probit` generate data from the normal model and the probit
  regression model, respectively.

## Example

We create a population of size $N = 50,000$ where our random variable of
interest $X$ and selection variable $Z$ jointly follow a bivariate
normal distribution with means $\mu_x = 10$, $\mu_z = 0$, standard
deviations $\sigma_x = 4$, $\sigma_z = 3$ and correlation $\rho = 0.5$.

For each population unit $l = 1, 2, \ldots, N$, the inclusion
probabilities are generated as $\pi_l = \Phi(b_0 + b_1 Z_l)$, where
$\Phi$ is the cumulative distribution function of a standard normal
distribution, $b_0 = -1.8$ and $b_1 = 0.1$.

Draw a sample of size $n = 500$ from the finite population with these
inclusion probabilities. Let
$w_{n,i} = n {\pi_i}^{-1} / \sum_{j=1}^n {\pi_j}^{-1}$ denote the scaled
sampling weight of sampled unit $i$, where $\pi_i$ denotes its inclusion
probability.

Our objective in this setting is to infer about the population mean,
$\mu_x$.

**Data generation** :

``` r
library(tidyverse)
source("data_generate.R")
source("normal_model.R")

set.seed(25)
# Specify the population and sample sizes
N = 50000; n = 500

# Specify the parameters to generate the data
mu_x = 10; mu_z = 0; s_x = 4; s_z = 3; rho = 0.5
b0 = -1.8; b1 =0.1;  M = 2000

# generate the data from the normal model
df = data_normal(mu_x = mu_x, mu_z = mu_z, s_x = s_x, s_z = s_z, rho = rho,
                 b0 = b0, b1 = b1, n = n, N = N)
```

**Get the point and interval estimates from the S-WLB method** :

``` r
est.SWLB = SWLB.normal(y = df$x, w = df$w, B = M)
est.SWLB
```

    ## $mu.SWLB
    ## [1] 9.946979
    ## 
    ## $sigma.sq.SWLB
    ## [1] 0.06501214
    ## 
    ## $CI.SWLB
    ##     2.5%    97.5% 
    ##  9.41759 10.40478

**Demonstrate the performance in point estimation and coverage
probability of uncertainty intervals** :

We perform $100$ Monte Carlo simulations to compare the performance in
point estimation and coverage probability of uncertainty intervals
across the different methods.

``` r
# Specify the number of Monte Carlo simulations
nRep = 100
# Specify the methods for comparison
methods = c("PMLE", "SWLB", "BPPE", "UBE")

# Matrices to store the point estimates and coverage indicators
est = CovI = matrix(NA, nrow = nRep, ncol = length(methods))
colnames(est) = colnames(CovI) = methods

set.seed(25)
for(r in 1:nRep){
  
  # generate the data
  df = data_normal(mu_x = mu_x, mu_z = mu_z, s_x = s_x, s_z = s_z, rho = rho,
                   b0 = b0, b1 = b1, n = n, N = N)
  
  # PMLE
  est.PMLE = PMLE.normal(y = df$x, w = df$w)
  est[r, "PMLE"] = est.PMLE$mu.PMLE
  CovI[r, "PMLE"] = 1*(mu_x >= est.PMLE$CI.PMLE[1] && mu_x <= est.PMLE$CI.PMLE[2])
  
  # UBE
  est.UBE = UBE.normal(y = df$x)
  est[r, "UBE"] = est.UBE$mu.UBE
  CovI[r, "UBE"] = 1*(mu_x >= est.UBE$CI.UBE[1] && mu_x <= est.UBE$CI.UBE[2])
  
  # BPPE
  est.BPPE = BPPE.normal(y = df$x, w = df$w)
  est[r, "BPPE"] = est.BPPE$mu.BPPE
  CovI[r, "BPPE"] = 1*(mu_x >= est.BPPE$CI.BPPE[1] && mu_x <= est.BPPE$CI.BPPE[2])
  
  # SWLB
  est.SWLB = SWLB.normal(y = df$x, w = df$w, B = M)
  est[r, "SWLB"] = est.SWLB$mu.SWLB
  CovI[r, "SWLB"] = 1*(mu_x >= est.SWLB$CI.SWLB[1] && mu_x <= est.SWLB$CI.SWLB[2])
}
```

Plot showing the point estimates:

``` r
# all plots follow a pre-specified theme stored in "themegg", 
# which can be found in the README.Rmd file.

# data frame for plotting the point estimates
df.est = data.frame(est = c(est), method = rep(methods, each = nRep))
df.est$method = factor(df.est$method, levels = methods)

# Boxplots of the point estimates 
g.est = ggplot(df.est, aes(x = method, y = est, color = method)) +
  geom_boxplot(outlier.size = 0.2) + 
  geom_hline(yintercept = mu_x, linetype = "dotted", color = "darkblue") +
  scale_color_manual(values = c( "brown2", "forestgreen","orange", "grey30")) +
  xlab("Methods") + ylab("Estimate of the mean") + themegg
g.est
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Plot showing the coverage probability of the $95$% intervals:

``` r
# all plots follow a pre-specified theme stored in "themegg", 
# which can be found in the README.Rmd file.

# Calculate the coverage probability
CovP = colMeans(CovI)

# data frame for plotting the coverage probability
df.cov = data.frame(covP = CovP, method = methods)
df.cov$method = factor(df.cov$method, levels = methods)

# Barplots of the coverage probability
g.cov = ggplot(df.cov, aes(x = method, y = covP, color = method, fill = method)) +
  geom_bar(stat="identity", width = 0.5, alpha = 0.6)+
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "darkblue") +
  scale_color_manual(values = c( "brown2", "forestgreen","orange", "grey30")) +
  scale_fill_manual(values = c( "brown2", "forestgreen","orange", "grey30")) +
  xlab("Methods") + ylab("Coverage probability") + themegg
g.cov
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
