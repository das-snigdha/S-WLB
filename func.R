library(mvtnorm)
library(foreach)
library(doParallel)
library(doRNG)

data_normal <- function(mu_x, mu_y, s_x, s_y, rho,
                        beta0, beta1, N,
                        pi_min = 0.02, pi_max = 0.98,
                        pilot.survey = FALSE, c_pi = 0.2) {
  
  mu <- c(mu_y, mu_x)
  
  S <- matrix(
    c(s_y^2, rho * s_x * s_y,
      rho * s_x * s_y, s_x^2),
    nrow = 2
  )
  
  data.pop <- mvtnorm::rmvnorm(n = N, mean = mu, sigma = S)
  colnames(data.pop) <- c("Y", "X")
  
  sel.prob <- pnorm(beta0 + beta1 * data.pop[, "X"])
  # sel.prob <- pmin(pmax(sel.prob, pi_min), pi_max)
  
  sel.ind <- rbinom(N, size = 1, prob = sel.prob)
  y <- data.pop[sel.ind == 1, "Y"]
  n_samp = length(y)
  
  raw.wts <- 1 / sel.prob[sel.ind == 1]
  sample.wts <- n_samp * raw.wts / sum(raw.wts)
  
  df.pop <- data.frame(
    Y = data.pop[, "Y"],
    X = data.pop[, "X"],
    sel.prob = sel.prob
  )
  
  sel.iid.ind = sample(1:N, size = n_samp, replace = FALSE)
  y.iid = data.pop[sel.iid.ind, "Y"]
  
  if(pilot.survey == TRUE){
    # sel.prob.pilot = pnorm(beta00 + beta11 * data.pop[, "X"])
    
    sel.prob.pilot = c_pi * sel.prob
    sel.ind.pilot <- rbinom(N, size = 1, prob = sel.prob.pilot)
    y.pilot <- data.pop[sel.ind.pilot == 1, "Y"]
    n.pilot = length(y.pilot)
    
    raw.wts.pilot <- 1 / sel.prob.pilot[sel.ind.pilot == 1]
    sample.wts.pilot <- n.pilot * raw.wts.pilot / sum(raw.wts.pilot)
    
    return(list(
      y = y,
      y.iid = y.iid,
      w = sample.wts,
      y.pilot = y.pilot,
      w.pilot = sample.wts.pilot,
      n.pilot = n.pilot,
      n_samp = n_samp,
      population = df.pop
    ))
  } else{
    return(list(
      y = y,
      y.iid = y.iid,
      w = sample.wts,
      n_samp = n_samp,
      population = df.pop
    ))
  } 
}


WFPBB = function(y, w, N){
  
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


oracle_iid_bayes <- function(y, s_y = 1, m0 = 0.35, tau = 0.2, B = 5000) {
  n <- length(y)
  V <- 1 / (n / s_y^2 + 1 / tau^2)
  m <- V * (sum(y) / s_y^2 + m0 / tau^2)
  
  draws <- rnorm(B, mean = m, sd = sqrt(V))
  list(mean = m, var = V, draws = draws)
}

pmle_mean <- function(y, w_tilde) {
  sum(w_tilde * y) / sum(w_tilde)
}

swlb_draws <- function(y, w_tilde, B = 5000) {
  n <- length(y)
  draws <- numeric(B)
  
  for (b in seq_len(B)) {
    G <- rexp(n, rate = 1)
    draws[b] <- sum(G * w_tilde * y) / sum(G * w_tilde)
  }
  
  draws
}

spb_draws <- function(y, w_tilde, s_y = 1, m0 = 0.35, tau = 0.2,
                      w0 = 1, B = 5000) {
  n <- length(y)
  draws <- numeric(B)
  
  for (b in seq_len(B)) {
    G <- rexp(n, rate = 1)
    
    num <- sum(G * w_tilde * y) / s_y^2 + w0 * m0 / tau^2
    den <- sum(G * w_tilde) / s_y^2 + w0 / tau^2
    
    draws[b] <- num / den
  }
  
  draws
}

bppe_draws <- function(y, w_tilde, s_y = 1, m0 = 0.35, tau = 0.2,
                       B = 5000) {
  V <- 1 / (sum(w_tilde) / s_y^2 + 1 / tau^2)
  m <- V * (sum(w_tilde * y) / s_y^2 + m0 / tau^2)
  
  rnorm(B, mean = m, sd = sqrt(V))
}

wfpbb_draws <- function(y, w_tilde, N,
                        s_y = 1, m0 = 0.35, tau = 0.2,
                        B = 5000) {
  
  draws <- numeric(B)
  w_raw = N*w_tilde/sum(w_tilde)
  
  for (b in seq_len(B)) {
    
    z <- WFPBB(y = y, w = w_raw, N = N)
    n_z <- length(z)
    
    V <- 1 / (n_z / s_y^2 + 1 / tau^2)
    m <- V * (sum(z) / s_y^2 + m0 / tau^2)
    
    draws[b] <- rnorm(1, mean = m, sd = sqrt(V))
  }
  
  draws
}

wasserstein_1d_normal_emp <- function(draws, oracle_draws) {
  mean((sort(draws) - sort(oracle_draws))^2)
}

estimate_w0_pompe <- function(y, w_tilde, s_y = 1) {
  mu_hat <- sum(w_tilde * y) / sum(w_tilde)
  Iw_hat <- mean(w_tilde^2 * ((y - mu_hat) / s_y^2)^2)
  I_hat <- 1 / s_y^2
  # sqrt(Iw_hat / I_hat)
  Iw_hat / I_hat
}

ci_draws <- function(draws, level = 0.95) {
  alpha <- 1 - level
  as.numeric(quantile(draws, probs = c(alpha / 2, 1 - alpha / 2)))
}

ise_to_oracle <- function(draws, oracle_mean, oracle_var, grid = NULL) {
  if (is.null(grid)) {
    s <- sqrt(oracle_var)
    grid <- seq(oracle_mean - 6 * s, oracle_mean + 6 * s, length.out = 1000)
  }
  
  oracle_dens <- dnorm(grid, mean = oracle_mean, sd = sqrt(oracle_var))
  dens <- density(draws, from = min(grid), to = max(grid), n = length(grid))
  dx <- diff(grid)[1]
  
  mean((dens$y - oracle_dens)^2) 
}

ise_normal_to_oracle <- function(mean_method, var_method, oracle_mean, oracle_var, grid = NULL) {
  if (is.null(grid)) {
    s <- max(sqrt(oracle_var), sqrt(var_method))
    center <- oracle_mean
    grid <- seq(center - 6 * s, center + 6 * s, length.out = 1000)
  }
  
  oracle_dens <- dnorm(grid, mean = oracle_mean, sd = sqrt(oracle_var))
  method_dens <- dnorm(grid, mean = mean_method, sd = sqrt(var_method))
  dx <- diff(grid)[1]
  
  mean((method_dens - oracle_dens)^2)
}


run_one_rep <- function(
    N = 10000,
    mu_x = 0, mu_y = 0, s_x = 3, s_y = 4,
    beta0 = -2.7, beta1 = 0.1, rho = 0.8,
    m0 = 0.25, tau = 0.2,
    B = 5000, est.w0 = TRUE,
    include_bwe = TRUE,
    pilot.survey = FALSE, c_pi = NULL
){
  
  df <- data_normal(
    mu_x = mu_x, mu_y = mu_y,
    s_x = s_x, s_y = s_y, rho = rho,
    beta0 = beta0, beta1 = beta1,
    N = N,
    pilot.survey = pilot.survey,
    c_pi = c_pi
  )
  
  y     <- df$y
  y.iid <- df$y.iid
  wt    <- df$w
  
  if (pilot.survey == TRUE) {
    pilot <- pilot_prior(
      y_p = df$y.pilot,
      w_p = df$w.pilot,
      m00 = m0,
      tau00 = tau,
      s_y = s_y
    )
    
    m_p <- pilot$m_p
    tau_p <- sqrt(pilot$V_p)
    n_pilot_out <- length(df$y.pilot)
  } else {
    m_p <- m0
    tau_p <- tau
    n_pilot_out <- NA
  }
  
  ## Oracle iid Bayes
  oracle <- oracle_iid_bayes(
    y.iid,
    s_y = s_y,
    m0 = m_p,
    tau = tau_p,
    B = B
  )
  
  oracle_mean <- oracle$mean
  oracle_var  <- oracle$var
  ci_oracle   <- ci_draws(oracle$draws)
  
  ## PMLE
  pmle <- pmle_mean(y, wt)
  
  Iw_hat <- mean(wt^2 * ((y - pmle) / s_y^2)^2)
  pmle_var <- s_y^4 * Iw_hat / length(y)
  
  ci_pmle <- c(
    pmle - qnorm(0.975) * sqrt(pmle_var),
    pmle + qnorm(0.975) * sqrt(pmle_var)
  )
  
  pmle_draws <- rnorm(B, mean = pmle, sd = sqrt(pmle_var))
  
  ## S-WLB
  swlb <- swlb_draws(y, wt, B)
  ci_swlb <- ci_draws(swlb)
  
  ## S-PB
  if (est.w0 == TRUE) {
    w0_hat <- estimate_w0_pompe(y, wt, s_y)
  } else {
    w0_hat <- 1
  }
  
  spb <- spb_draws(
    y, wt,
    s_y = s_y,
    m0 = m_p,
    tau = tau_p,
    w0 = w0_hat,
    B = B
  )
  
  ci_spb <- ci_draws(spb)
  
  ## BPPE
  bppe <- bppe_draws(
    y, wt,
    s_y = s_y,
    m0 = m_p,
    tau = tau_p,
    B = B
  )
  
  ci_bppe <- ci_draws(bppe)
  bppe_var <- 1 / (sum(wt) / s_y^2 + 1 / tau_p^2)
  
  ## WFPBB / BWE (optional)
  if (include_bwe) {
    
    wfpbb <- wfpbb_draws(
      y = y,
      w_tilde = wt,
      N = N,
      s_y = s_y,
      m0 = m_p,
      tau = tau_p,
      B = B
    )
    
    ci_wfpbb <- ci_draws(wfpbb)
    
  } else {
    
    wfpbb <- rep(NA_real_, B)
    ci_wfpbb <- c(NA_real_, NA_real_)
    
  }
  
  ## UBE
  n.samp <- length(y)
  ube <- bppe_draws(
    y, rep(1, n.samp),
    s_y = s_y,
    m0 = m_p,
    tau = tau_p,
    B = B
  )
  
  ci_ube <- ci_draws(ube)
  ube_var <- 1 / (n.samp / s_y^2 + 1 / tau_p^2)
  
  ## Prior-induced shifts
  iid_mean <- mean(y.iid)
  
  delta_oracle <- oracle_mean - iid_mean
  delta_pmle   <- 0
  delta_swlb   <- mean(swlb) - pmle
  delta_spb    <- mean(spb) - pmle
  delta_bppe   <- mean(bppe) - pmle
  delta_wfpbb  <- mean(wfpbb) - pmle
  delta_ube    <- mean(ube) - pmle
  
  summary_row <- data.frame(
    
    n_survey = length(y),
    w0_hat = w0_hat,
    
    n.pilot = n_pilot_out,
    prior_mean = m_p,
    prior_sd = tau_p,
    
    oracle_mean = oracle_mean,
    oracle_var = oracle_var,
    oracle_bias = oracle_mean - mu_y,
    oracle_cover = as.integer(
      ci_oracle[1] <= mu_y && mu_y <= ci_oracle[2]
    ),
    oracle_length = diff(ci_oracle),
    
    pmle = pmle,
    swlb_mean = mean(swlb),
    spb_mean = mean(spb),
    bppe_mean = mean(bppe),
    wfpbb_mean = mean(wfpbb),
    ube_mean = mean(ube),
    
    bias_pmle = pmle - mu_y,
    bias_swlb = mean(swlb) - mu_y,
    bias_spb = mean(spb) - mu_y,
    bias_bppe = mean(bppe) - mu_y,
    bias_wfpbb = mean(wfpbb) - mu_y,
    bias_ube = mean(ube) - mu_y,
    
    cover_pmle = as.integer(
      ci_pmle[1] <= mu_y && mu_y <= ci_pmle[2]
    ),
    cover_swlb = as.integer(
      ci_swlb[1] <= mu_y && mu_y <= ci_swlb[2]
    ),
    cover_spb = as.integer(
      ci_spb[1] <= mu_y && mu_y <= ci_spb[2]
    ),
    cover_bppe = as.integer(
      ci_bppe[1] <= mu_y && mu_y <= ci_bppe[2]
    ),
    cover_wfpbb = as.integer(
      ci_wfpbb[1] <= mu_y && mu_y <= ci_wfpbb[2]
    ),
    cover_ube = as.integer(
      ci_ube[1] <= mu_y && mu_y <= ci_ube[2]
    ),
    
    length_pmle = diff(ci_pmle),
    length_swlb = diff(ci_swlb),
    length_spb = diff(ci_spb),
    length_bppe = diff(ci_bppe),
    length_wfpbb = diff(ci_wfpbb),
    length_ube = diff(ci_ube),
    
    dist_pmle = abs(mean(pmle_draws) - oracle_mean),
    dist_swlb = abs(mean(swlb) - oracle_mean),
    dist_spb = abs(mean(spb) - oracle_mean),
    dist_bppe = abs(mean(bppe) - oracle_mean),
    dist_wfpbb = abs(mean(wfpbb) - oracle_mean),
    dist_ube = abs(mean(ube) - oracle_mean),
    
    w2_swlb = wasserstein_1d_normal_emp(
      swlb, oracle$draws
    ),
    w2_spb = wasserstein_1d_normal_emp(
      spb, oracle$draws
    ),
    w2_bppe = wasserstein_1d_normal_emp(
      bppe, oracle$draws
    ),
    w2_pmle = wasserstein_1d_normal_emp(
      pmle_draws, oracle$draws
    ),
    w2_wfpbb = if (include_bwe) {
      wasserstein_1d_normal_emp(wfpbb, oracle$draws)
    } else {
      NA_real_
    },
    w2_ube = wasserstein_1d_normal_emp(
      ube, oracle$draws
    ),
    
    ise_pmle = ise_normal_to_oracle(
      mean(pmle_draws), pmle_var,
      oracle_mean, oracle_var
    ),
    ise_swlb = ise_to_oracle(
      swlb, oracle_mean, oracle_var
    ),
    ise_spb = ise_to_oracle(
      spb, oracle_mean, oracle_var
    ),
    ise_bppe = ise_normal_to_oracle(
      mean(bppe), bppe_var,
      oracle_mean, oracle_var
    ),
    ise_wfpbb = if (include_bwe) {
      ise_to_oracle(wfpbb, oracle_mean, oracle_var)
    } else {
      NA_real_
    },
    ise_ube = ise_normal_to_oracle(
      mean(ube), ube_var,
      oracle_mean, oracle_var
    ),
    
    delta_oracle = delta_oracle,
    delta_pmle = delta_pmle,
    delta_swlb = delta_swlb,
    delta_spb = delta_spb,
    delta_bppe = delta_bppe,
    delta_wfpbb = delta_wfpbb,
    delta_ube = delta_ube,
    
    shift_dist_pmle = abs(delta_pmle - delta_oracle),
    shift_dist_swlb = abs(delta_swlb - delta_oracle),
    shift_dist_spb = abs(delta_spb - delta_oracle),
    shift_dist_bppe = abs(delta_bppe - delta_oracle),
    shift_dist_wfpbb = abs(delta_wfpbb - delta_oracle),
    shift_dist_ube = abs(delta_ube - delta_oracle)
  )
  
  draws <- list(
    oracle = oracle$draws,
    pmle = pmle_draws,
    swlb = swlb,
    spb = spb,
    bppe = bppe,
    wfpbb = wfpbb,
    ube = ube
  )
  
  list(summary = summary_row, draws = draws)
}

run_simulation_parallel <- function(R = 200, seed = 123,
                                    n_cores = 50,
                                    print.every = 20, ...) {
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  })
  
  out <- foreach(
    r = seq_len(R),
    .options.RNG = seed,
    .packages = c("mvtnorm"),
    .export = c(
      "run_one_rep", "data_normal", "pilot_prior",
      "WFPBB",
      "oracle_iid_bayes", "pmle_mean",
      "swlb_draws", "spb_draws", "bppe_draws", "wfpbb_draws",
      "estimate_w0_pompe", "ci_draws",
      "wasserstein_1d_normal_emp",
      "ise_to_oracle", "ise_normal_to_oracle"
    )
  ) %dorng% {
    
    tmp <- run_one_rep(...)
    
    if (r %% print.every == 0) {
      cat("Finished", r, "replicates\n")
    }
    
    tmp
  }
  
  summaries <- lapply(out, function(x) x$summary)
  draws_list <- lapply(out, function(x) x$draws)
  
  list(
    res = do.call(rbind, summaries),
    draws_list = draws_list
  )
}


pilot_prior <- function(
    y_p, w_p, m00 = 0.25, tau00 = 40, s_y = 1
) {
  
  V_post <- 1 / (sum(w_p) / s_y^2 + 1 / tau00^2)
  
  m_p <- V_post * (
    sum(w_p * y_p) / s_y^2 +
      m00 / tau00^2
  )
  
  mu_p <- sum(w_p * y_p) / sum(w_p)
  
  Iw_hat <- mean(w_p^2 * ((y_p - m_p) / s_y^2)^2)
  V_p <- s_y^4 * Iw_hat / length(y_p)
  
  list(
    m_p = m_p,
    V_p = V_p,
    n_pilot = length(y_p),
    y_pilot = y_p,
    w_pilot = w_p
  )
}