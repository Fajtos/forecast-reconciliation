rm(list = ls()); gc()

library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tidyr)
has_crps <- requireNamespace("scoringRules", quietly = TRUE)
if (has_crps) library(scoringRules)

# DATA GENERATION & METRICS 

simulate_forecast_matrix <- function(T = 100, H = 5, phi = 0.8, sigma_eta = 0.3, sigma_y = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu <- numeric(T + H)
  mu[1] <- rnorm(1, 0, sigma_eta / sqrt(1 - phi^2))
  for (t in 2:(T + H)) mu[t] <- phi * mu[t - 1] + rnorm(1, 0, sigma_eta)
  mu <- mu / sd(mu) * 0.5
  
  delta <- rnorm(H, 0, 0.1); delta <- delta - mean(delta)
  tau   <- 0.05 + 0.02 * (1:H) + 0.005 * sin(pi * (1:H) / H)
  
  forecast_matrix <- matrix(NA_real_, nrow = T, ncol = H)
  for (t in 1:T) {
    for (h in 1:min(H, T + H - t)) {
      eps <- rnorm(1, 0, tau[h])
      if (h > 1) eps <- 0.7 * eps + 0.3 * rnorm(1, 0, tau[h])
      if (runif(1) < 0.05) eps <- eps * 3
      forecast_matrix[t, h] <- mu[t + h] + delta[h] + eps
    }
  }
  y <- mu[1:T] + rnorm(T, 0, sigma_y)
  list(T = T, H = H, mu = mu[1:T], forecasts = forecast_matrix, y = y, true_bias = delta, true_tau = tau)
}

simulate_forecast_matrix_tar <- function(T = 100, H = 5, phi = 0.8, sigma_eta = 0.3, sigma_y = 0.15, tar_strength = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu <- numeric(T + H)
  mu[1] <- rnorm(1, 0, sigma_eta / sqrt(1 - phi^2))
  for (t in 2:(T + H)) mu[t] <- phi * mu[t - 1] + rnorm(1, 0, sigma_eta)
  mu <- mu / sd(mu) * 0.5
  
  delta <- rnorm(H, 0, 0.1); delta <- delta - mean(delta)
  tau   <- 0.05 + 0.02 * (1:H) + 0.005 * sin(pi * (1:H) / H)
  
  forecast_matrix <- matrix(NA_real_, nrow = T, ncol = H)
  for (t in 1:T) {
    for (h in 1:min(H, T + H - t)) {
      prev <- if (t + h - 1 >= 1) mu[t + h - 1] else 0
      thr  <- if (prev > 0) 1 else -1
      eps  <- rt(1, df = 5) * tau[h] + tar_strength * thr * tau[h]
      forecast_matrix[t, h] <- mu[t + h] + delta[h] + eps
    }
  }
  y <- mu[1:T] + rnorm(T, 0, sigma_y)
  list(T = T, H = H, mu = mu[1:T], forecasts = forecast_matrix, y = y, true_bias = delta, true_tau = tau)
}

#  metric matching target = origin + horizon
calculate_incoherence <- function(forecast_matrix) {
  T <- nrow(forecast_matrix); H <- ncol(forecast_matrix)
  total_sse <- 0; count <- 0
  for (s in seq_len(T)) {
    vals <- numeric(0)
    if (s > 1) { # To avoid R sequence bug 1:0
      for (i in max(1, s - H):(s - 1)) {
        h <- s - i
        if (h >= 1 && h <= H && !is.na(forecast_matrix[i, h])) {
          vals <- c(vals, forecast_matrix[i, h])
        }
      }
    }
    if (length(vals) > 1) {
      for (k in 1:(length(vals) - 1)) {
        for (j in (k + 1):length(vals)) {
          total_sse <- total_sse + (vals[k] - vals[j])^2
          count <- count + 1
        }
      }
    }
  }
  if (count == 0) return(0)
  total_sse / count
}

# STAN PREPARATION 

if (!dir.exists("results")) dir.create("results")

mod <- cmdstan_model("reconcile.stan", force_recompile = TRUE)

prepare_stan_data <- function(fm, y_vec = NULL, y_idx = NULL, prior_c_mean = 0, prior_c_sd = 0.5,
                              prior_phi_mean = 0.8, prior_phi_sd = 0.1, prior_sigma_eta_shape = 3, prior_sigma_eta_rate = 2,
                              prior_delta_sd = 0.5, prior_tau_shape = 5, prior_tau_rate = 2,
                              prior_sigma_y_shape = 2, prior_sigma_y_rate = 0.1) {
  T <- nrow(fm); H <- ncol(fm)
  origin_idx <- integer(0); horizon_idx <- integer(0); forecast <- numeric(0)
  for (t in seq_len(T)) for (h in seq_len(H)) if (!is.na(fm[t, h])) {
    origin_idx <- c(origin_idx, t); horizon_idx <- c(horizon_idx, h); forecast <- c(forecast, fm[t, h])
  }
  list(
    T = T, H = H, N = length(forecast), origin_idx = origin_idx, horizon_idx = horizon_idx, forecast = forecast,
    prior_phi_mean = prior_phi_mean, prior_phi_sd = prior_phi_sd, prior_sigma_eta_shape = prior_sigma_eta_shape,
    prior_sigma_eta_rate = prior_sigma_eta_rate, prior_delta_sd = prior_delta_sd, prior_tau_shape = prior_tau_shape,
    prior_tau_rate = prior_tau_rate, prior_c_mean = prior_c_mean, prior_c_sd = prior_c_sd,
    T_y = if (is.null(y_vec)) 0L else length(y_vec), y_obs = if (is.null(y_vec)) numeric(0) else as.numeric(y_vec),
    y_idx = if (is.null(y_idx)) integer(0) else as.integer(y_idx),
    prior_sigma_y_shape = prior_sigma_y_shape, prior_sigma_y_rate = prior_sigma_y_rate
  )
}

# RUN ONE REPLICATION

run_one <- function(rep_id, T = 100, H = 5, stress = FALSE, tar_strength = 2.0,
                    chains = 2, iter_warmup = 500, iter_sampling = 500) {
  sim <- if (stress) simulate_forecast_matrix_tar(T, H, tar_strength = tar_strength, seed = rep_id)
  else simulate_forecast_matrix(T, H, seed = rep_id)
  
  sd_stan <- prepare_stan_data(sim$forecasts, y_vec = sim$y, y_idx = seq_len(T), prior_c_sd = 0.5)
  
  fit <- mod$sample(data = sd_stan, chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                    parallel_chains = chains, refresh = 0, adapt_delta = 0.95, max_treedepth = 12)
  
  mu_mean  <- colMeans(fit$draws("mu", format = "matrix"))
  est_bias <- colMeans(fit$draws("delta", format = "matrix"))
  
  rec_debiased <- matrix(NA_real_, T, H)
  rec_coherent <- matrix(NA_real_, T, H)
  for (t in 1:T) for (h in 1:H) if (!is.na(sim$forecasts[t, h])) {
    rec_debiased[t, h] <- mu_mean[t + h] + est_bias[h]
    rec_coherent[t, h] <- mu_mean[t + h]
  }
  
  raw_inc      <- calculate_incoherence(sim$forecasts)
  inc_debiased <- calculate_incoherence(rec_debiased)
  inc_coherent <- calculate_incoherence(rec_coherent)
  
  red_debiased <- if (raw_inc > 0) (raw_inc - inc_debiased) / raw_inc * 100 else 0
  red_coherent <- if (raw_inc > 0) (raw_inc - inc_coherent) / raw_inc * 100 else 0
  
  cor_mu <- cor(mu_mean[1:T], sim$mu)
  
  y_rep <- fit$draws("y_rep", format = "matrix")
  lo <- apply(y_rep, 2, quantile, 0.05); hi <- apply(y_rep, 2, quantile, 0.95)
  cov90 <- mean(sim$y >= lo & sim$y <= hi)
  crps  <- if (has_crps) mean(crps_sample(sim$y, t(y_rep))) else NA_real_
  
  params <- c("c", "phi", "sigma_eta_sq", "delta_raw", "tau_sq", "sigma_y_sq")
  summ   <- fit$summary(variables = params)
  conv   <- all(summ$rhat < 1.05, na.rm = TRUE)
  
  list(rep = rep_id, raw_inc = raw_inc, inc_debiased = inc_debiased, inc_coherent = inc_coherent,
       reduction_debiased = red_debiased, reduction_coherent = red_coherent, mu_corr = cor_mu,
       c_est = mean(fit$draws("c", format = "matrix")), coverage90 = cov90, crps = crps, converged = conv)
}

# BATCH RUNNER 

run_batch <- function(R, T = 100, H = 5, stress = FALSE, tar_strength = 2.0, tag = "baseline") {
  results <- vector("list", R)
  for (r in 1:R) {
    cat(sprintf("[%s] rep %d/%d\n", tag, r, R))
    results[[r]] <- tryCatch(
      run_one(r, T, H, stress = stress, tar_strength = tar_strength),
      error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL }
    )
  }
  Filter(Negate(is.null), results)
}

summarise_results <- function(results, tag) {
  if (length(results) == 0) stop(sprintf("[%s] No results to summarise.", tag))
  to_num <- function(x) as.numeric(unlist(x))
  
  df <- data.frame(
    rep = to_num(sapply(results, `[[`, "rep")),
    reduction_debiased = to_num(sapply(results, `[[`, "reduction_debiased")),
    reduction_coherent = to_num(sapply(results, `[[`, "reduction_coherent")),
    mu_corr = to_num(sapply(results, `[[`, "mu_corr")),
    c_est = to_num(sapply(results, `[[`, "c_est")),
    coverage90 = to_num(sapply(results, `[[`, "coverage90")),
    converged = to_num(sapply(results, `[[`, "converged"))
  )
  
  out_file <- sprintf("results/simulation_summary_%s.csv", tag)
  tryCatch({
    write.csv(df, out_file, row.names = FALSE)
  }, error = function(e) {
    warning(sprintf("Could not write %s (Is it open in Excel?). Results are still in memory.", out_file))
  })
  
  cat(sprintf("\n[%s] N = %d\n", tag, nrow(df)))
  cat(sprintf("[%s] reduction (de-biased): %.1f%%\n", tag, mean(df$reduction_debiased)))
  cat(sprintf("[%s] reduction (coherent):  %.1f%%\n", tag, mean(df$reduction_coherent)))
  cat(sprintf("[%s] mu corr:                %.3f\n", tag, mean(df$mu_corr)))
  cat(sprintf("[%s] coverage (nom. 90%%):   %.1f%%\n", tag, 100 * mean(df$coverage90)))
  cat(sprintf("[%s] convergence:            %.0f%%\n", tag, 100 * mean(df$converged)))
  return(df)
}

#  EXECUTION 

cat("\n SMOKE TEST (5 reps) \n")
res_quick <- run_batch(5, tag = "quick")
df_quick  <- summarise_results(res_quick, "quick")

cat("\n BASELINE (100 reps) \n")
res_main <- run_batch(100, tag = "baseline")
df_main  <- summarise_results(res_main, "baseline")

cat("\n STRESS TEST (50 reps) \n")
res_stress <- run_batch(50, tag = "stress", stress = TRUE)
df_stress  <- summarise_results(res_stress, "stress")

# PLOTTING 

p1 <- ggplot(df_main, aes(reduction_debiased)) + geom_histogram(bins=20, fill="steelblue", colour="black", alpha=.8) + geom_vline(xintercept=mean(df_main$reduction_debiased), linetype="dashed", colour="red") + labs(x="Incoherence reduction (%) [de-biased]", y="Frequency") + theme_minimal()
ggsave("results/hist_red_debiased.png", p1, width=7, height=5, dpi=150)

p2 <- ggplot(df_main, aes(reduction_coherent)) + geom_histogram(bins=20, fill="darkorange", colour="black", alpha=.8) + labs(x="Incoherence reduction (%) [coherent]", y="Frequency") + theme_minimal()
ggsave("results/hist_red_coherent.png", p2, width=7, height=5, dpi=150)

p3 <- ggplot(df_main, aes(mu_corr)) + geom_histogram(bins=20, fill="darkgreen", colour="black", alpha=.8) + geom_vline(xintercept=mean(df_main$mu_corr), linetype="dashed", colour="red") + labs(x=expression("Correlation with true "*mu), y="Frequency") + theme_minimal()
ggsave("results/hist_mu_corr.png", p3, width=7, height=5, dpi=150)

p4 <- ggplot(df_main, aes(reduction_debiased, mu_corr)) + geom_point(alpha=.7, colour="darkgreen") + geom_smooth(method="lm", se=FALSE, colour="red", linetype="dashed") + labs(x="Reduction (%) [de-biased]", y=expression("Correlation with true "*mu)) + theme_minimal()
ggsave("results/scatter_red_vs_corr.png", p4, width=7, height=5, dpi=150)

p5 <- ggplot(df_main, aes(coverage90)) + geom_histogram(bins=20, fill="purple", colour="black", alpha=.8) + geom_vline(xintercept=0.90, linetype="dashed", colour="red") + labs(x="Coverage of nominal 90% PI", y="Frequency") + theme_minimal()
ggsave("results/hist_coverage.png", p5, width=7, height=5, dpi=150)

p6 <- ggplot(df_stress, aes(mu_corr)) + geom_histogram(bins=20, fill="firebrick", colour="black", alpha=.8) + geom_vline(xintercept=mean(df_stress$mu_corr), linetype="dashed", colour="red") + labs(x=expression("Correlation with true "*mu*" [stress test]"), y="Frequency") + theme_minimal()
ggsave("results/hist_mu_corr_stress.png", p6, width=7, height=5, dpi=150)


# FINAL METRICS SUMMARY (FOR REPRODUCIBILITY)
cat("\n FINAL SUMMARY METRICS \n")

# Load the saved CSVs
base_df   <- read.csv("results/simulation_summary_baseline.csv")
stress_df <- read.csv("results/simulation_summary_stress.csv")

# Helper function for formatting
print_metric <- function(label, value, is_pct = FALSE, digits = 3) {
  if (is_pct) {
    cat(sprintf("%-30s : %5.1f%%\n", label, value * 100))
  } else {
    cat(sprintf("%-30s : %5.*f\n", label, digits, value))
  }
}

cat("\n BASELINE (AR(1)) \n")
print_metric("N Replications", nrow(base_df), digits = 0)
# REMOVED is_pct = TRUE for reductions because they are already 0-100
print_metric("Mean Reduction (de-biased)", mean(base_df$reduction_debiased), is_pct = FALSE, digits = 1) 
print_metric("Mean Reduction (coherent)", mean(base_df$reduction_coherent), is_pct = FALSE, digits = 1)
print_metric("Mean Correlation with true mu", mean(base_df$mu_corr), digits = 3)
# KEPT is_pct = TRUE for coverage because it is 0-1
print_metric("Mean Coverage (90% PI)", mean(base_df$coverage90), is_pct = TRUE) 
print_metric("Convergence Rate", mean(base_df$converged), is_pct = TRUE)

cat("\n STRESS TEST (TAR + Student-t) \n")
print_metric("N Replications", nrow(stress_df), digits = 0)
# REMOVED is_pct = TRUE for reductions
print_metric("Mean Reduction (de-biased)", mean(stress_df$reduction_debiased), is_pct = FALSE, digits = 1)
print_metric("Mean Reduction (coherent)", mean(stress_df$reduction_coherent), is_pct = FALSE, digits = 1)
print_metric("Mean Correlation with true mu", mean(stress_df$mu_corr), digits = 3)
print_metric("Mean Coverage (90% PI)", mean(stress_df$coverage90), is_pct = TRUE)
print_metric("Convergence Rate", mean(stress_df$converged), is_pct = TRUE)