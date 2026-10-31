# COMPLETE EMPIRICAL PIPELINE - FULL EXPORTS & PLOTS
rm(list = ls()); gc()
library(cmdstanr)
library(forecast)
library(quantmod)
library(scoringRules)
library(dplyr)
library(ggplot2)
library(tidyr)

# HELPER FUNCTIONS
prepare_stan_data <- function(fm, y_vec = NULL, y_idx = NULL, prior_c_mean = 0, prior_c_sd = 1,
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

calculate_incoherence <- function(forecast_matrix) {
  T <- nrow(forecast_matrix); H <- ncol(forecast_matrix)
  total_sse <- 0; count <- 0
  for (s in seq_len(T)) {
    vals <- numeric(0)
    if (s > 1) {
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

diag_recon <- function(fm) {
  T <- nrow(fm); H <- ncol(fm)
  rec <- matrix(NA, T, H)
  for (s in 1:T) {
    if (s > 1) {
      vals <- c(); idx <- list()
      start_i <- max(1, s - H); end_i <- s - 1
      if (start_i <= end_i) {
        for (i in start_i:end_i) {
          h <- s - i
          if (h >= 1 && h <= H && !is.na(fm[i, h])) { 
            vals <- c(vals, fm[i, h]); idx[[length(idx) + 1]] <- c(i, h) 
          }
        }
      }
      if (length(vals) > 0) {
        avg <- mean(vals)
        for (k in 1:length(idx)) rec[idx[[k]][1], idx[[k]][2]] <- avg
      }
    }
  }
  rec
}

# LOAD DATA

getSymbols("SPY", from = "2023-12-01", to = "2025-12-31", auto.assign = TRUE)
spy_ret <- na.omit(dailyReturn(SPY$SPY.Adjusted))
spy_vol <- na.omit(rollapply(spy_ret, 5, sd, align = "right") * sqrt(252))
y_all   <- as.numeric(spy_vol)
T_full  <- length(y_all)
T_test  <- 100
H       <- 5
y_oos   <- y_all[(T_full - T_test + 1):T_full]
vol_dates <- index(spy_vol)[(T_full - T_test + 1):T_full]

#GENERATE FORECASTS

fm_ar1 <- matrix(NA, T_test, H)
fm_ets <- matrix(NA, T_test, H)

for (t in 1:(T_test - H + 1)) {
  train_end <- T_full - T_test + t - 1
  train     <- y_all[1:train_end]
  fm_ar1[t, ] <- as.numeric(forecast(Arima(train, order = c(1,0,0)), h = H)$mean)
  fm_ets[t, ] <- as.numeric(forecast(ets(train), h = H)$mean)
}

if (!file.exists("lstm_forecasts.csv")) stop("lstm_forecasts.csv not found.")
fm_lstm <- as.matrix(read.csv("lstm_forecasts.csv"))
colnames(fm_lstm) <- NULL

#RECONCILIATION METHODS

mod <- cmdstan_model("reconcile.stan", force_recompile = TRUE)

run_bayesian <- function(fm) {
  s <- sd(fm, na.rm = TRUE); m <- mean(fm, na.rm = TRUE)
  fm_std <- (fm - m) / s; y_std  <- (y_oos - m) / s
  
  sd_stan <- prepare_stan_data(fm_std, y_vec = y_std, y_idx = seq_len(T_test), prior_c_sd = 1)
  fit <- mod$sample(data = sd_stan, chains = 4, iter_warmup = 1000, iter_sampling = 1000, 
                    parallel_chains = 4, refresh = 0, show_messages = FALSE)
  
  mu_mean <- colMeans(fit$draws("mu", format = "matrix")) * s + m
  delta   <- colMeans(fit$draws("delta", format = "matrix")) * s
  tau     <- colMeans(fit$draws("tau", format = "matrix")) * s
  y_rep   <- fit$draws("y_rep", format = "matrix") * s + m
  
  rec <- matrix(NA, T_test, H)
  for (t in 1:T_test) for (h in 1:H) if (!is.na(fm[t, h])) rec[t, h] <- mu_mean[t + h] + delta[h]
  list(rec = rec, y_rep = y_rep, mu = mu_mean, delta = delta, tau = tau, fit = fit)
}

cat("\nFitting Bayesian Models...\n")
res_ar1  <- run_bayesian(fm_ar1)
res_ets  <- run_bayesian(fm_ets)
res_lstm <- run_bayesian(fm_lstm)
diag_ar1 <- diag_recon(fm_ar1); diag_ets <- diag_recon(fm_ets); diag_lstm <- diag_recon(fm_lstm)

# METRICS & CSVS

calc_metrics <- function(y_true, y_rep, fm_raw, fm_rec_bayes, fm_rec_diag) {
  rmse_raw <- sqrt(mean((rowMeans(fm_raw, na.rm=T) - y_true)^2, na.rm=T))
  rmse_bayes <- sqrt(mean((rowMeans(fm_rec_bayes, na.rm=T) - y_true)^2, na.rm=T))
  rmse_diag  <- sqrt(mean((rowMeans(fm_rec_diag, na.rm=T) - y_true)^2, na.rm=T))
  
  crps_bayes <- mean(crps_sample(y_true, t(y_rep)))
  cov_bayes  <- mean(y_true >= apply(y_rep, 2, quantile, 0.05) & y_true <= apply(y_rep, 2, quantile, 0.95))
  
  inc_raw <- calculate_incoherence(fm_raw)
  inc_diag <- calculate_incoherence(fm_rec_diag)
  inc_bayes <- calculate_incoherence(fm_rec_bayes)
  
  data.frame(RMSE_Raw = rmse_raw, RMSE_Diag = rmse_diag, RMSE_Bayes = rmse_bayes, 
             CRPS_Bayes = crps_bayes, Coverage_Bayes = cov_bayes,
             Incoh_Raw = inc_raw, Incoh_Diag = inc_diag, Incoh_Bayes = inc_bayes)
}

results <- rbind(
  AR1  = calc_metrics(y_oos, res_ar1$y_rep, fm_ar1, res_ar1$rec, diag_ar1),
  ETS  = calc_metrics(y_oos, res_ets$y_rep, fm_ets, res_ets$rec, diag_ets),
  LSTM = calc_metrics(y_oos, res_lstm$y_rep, fm_lstm, res_lstm$rec, diag_lstm)
)
write.csv(results, "table_performance.csv")

# Horizon Diagnostics CSV
diag_ets_df <- data.frame(Model="ETS", Horizon=1:H, Bias=res_ets$delta, Reliability=res_ets$tau)
diag_lstm_df <- data.frame(Model="LSTM", Horizon=1:H, Bias=res_lstm$delta, Reliability=res_lstm$tau)
write.csv(rbind(diag_ets_df, diag_lstm_df), "table_diagnostics.csv", row.names=FALSE)

# VaR Backtest CSV
daily_vol_raw   <- rowMeans(fm_ar1, na.rm = T) / sqrt(252)
daily_vol_diag  <- rowMeans(diag_ar1, na.rm = T) / sqrt(252)
daily_vol_bayes <- rowMeans(res_ar1$rec, na.rm = T) / sqrt(252)
actual_ret      <- as.numeric(spy_ret[(length(spy_ret) - T_test + 1):length(spy_ret)])

var_results <- data.frame(
  Method = c("Expected (5%)", "Raw AR1", "Diagonal Baseline", "Bayesian Reconciled"),
  VaR_Breaches = c(5, 
                   sum(actual_ret < -1.645 * daily_vol_raw, na.rm = T),
                   sum(actual_ret < -1.645 * daily_vol_diag, na.rm = T),
                   sum(actual_ret < -1.645 * daily_vol_bayes, na.rm = T))
)
write.csv(var_results, "table_var.csv", row.names=FALSE)

# MCMC Convergence CSV
mcmc_diag <- res_lstm$fit$summary(variables = c("c", "phi", "sigma_eta_sq", "delta_raw", "tau_sq", "sigma_y_sq"))
write.csv(mcmc_diag[, c("variable", "rhat", "ess_bulk")], "table_mcmc.csv", row.names=FALSE)

cat("\nAll result CSVs written successfully.\n")

#PLOTTING

dir.create("figures", showWarnings = FALSE)

# Horizon Specific Diagnostics Plot
plot_diagnostics <- function(res, model_name) {
  df <- data.frame(Horizon = 1:H, Bias = res$delta, Reliability = res$tau)
  df_long <- pivot_longer(df, cols = c("Bias", "Reliability"), names_to = "Parameter", values_to = "Value")
  p <- ggplot(df_long, aes(x = Horizon, y = Value, fill = Parameter)) +
    geom_col(position = "dodge", alpha = 0.8) + facet_wrap(~ Parameter, scales = "free_y") +
    scale_x_continuous(breaks = 1:H) + labs(title = paste(model_name, "Horizon-Specific Parameters")) +
    theme_minimal() + theme(legend.position = "none")
  ggsave(paste0("figures/diagnostics_", model_name, ".png"), p, width = 8, height = 4, dpi = 150)
}
plot_diagnostics(res_ets, "ETS")
plot_diagnostics(res_lstm, "LSTM")

# Raw vs Reconciled Trajectories Plot
plot_trajectory <- function(fm_raw, res, model_name, start_idx = 70, end_idx = 100) {
  target_times <- start_idx:end_idx
  n_steps <- length(target_times)
  raw_vals <- numeric(n_steps); rec_vals <- numeric(n_steps)
  actual_vals <- y_oos[target_times]; dates_sub <- vol_dates[target_times]
  
  for (i in seq_along(target_times)) {
    t <- target_times[i]
    orig <- max(1, t - H + 1); h <- t - orig + 1
    raw_vals[i] <- fm_raw[orig, h]; rec_vals[i] <- res$rec[orig, h]
  }
  
  df <- data.frame(Date = dates_sub, Actual = actual_vals, Raw = raw_vals, Reconciled = rec_vals)
  p <- ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Actual, color = "Actual"), linewidth = 1) +
    geom_point(aes(y = Raw, color = "Raw Forecasts"), size = 2, alpha = 0.6) +
    geom_line(aes(y = Reconciled, color = "Reconciled"), linewidth = 1) +
    scale_color_manual(values = c("Actual" = "black", "Raw Forecasts" = "red", "Reconciled" = "blue")) +
    labs(title = paste(model_name, ": Raw vs Reconciled (Last 30 Days)"), x = "Date", y = "Volatility", color = "") + 
    theme_minimal()
  ggsave(paste0("figures/trajectory_", model_name, ".png"), p, width = 8, height = 4, dpi = 150)
}
plot_trajectory(fm_lstm, res_lstm, "LSTM")
plot_trajectory(fm_ets, res_ets, "ETS")

# Posterior Predictive Intervals Plot
plot_intervals <- function(res, model_name) {
  y_rep <- res$y_rep
  lo <- apply(y_rep, 2, quantile, 0.05); hi <- apply(y_rep, 2, quantile, 0.95); med <- apply(y_rep, 2, median)
  df <- data.frame(Date = vol_dates, Actual = y_oos, Median = med, Lower = lo, Upper = hi)
  p <- ggplot(df, aes(x = Date)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "90% PI"), alpha = 0.3) +
    geom_line(aes(y = Actual, color = "Actual"), linewidth = 0.8) +
    geom_line(aes(y = Median, color = "Posterior Median"), linewidth = 0.8, linetype = "dashed") +
    scale_fill_manual(values = c("90% PI" = "steelblue")) +
    scale_color_manual(values = c("Actual" = "black", "Posterior Median" = "red")) +
    labs(title = paste(model_name, ": Posterior Predictive Intervals"), x = "Date", y = "Volatility", fill = "", color = "") + 
    theme_minimal()
  ggsave(paste0("figures/intervals_", model_name, ".png"), p, width = 8, height = 4, dpi = 150)
}
plot_intervals(res_lstm, "LSTM")

#Distribution Check (Empirical vs Theoretical Normal) Plot
plot_distribution <- function(fm_raw, model_name) {
  # We use the raw forecasts to check the marginal distribution
  vals <- as.numeric(fm_raw)
  vals <- vals[!is.na(vals)]
  
  p <- ggplot(data.frame(x = vals), aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "grey70", colour = "black", alpha = 0.8) +
    stat_function(fun = dnorm, args = list(mean = mean(vals), sd = sd(vals)), 
                  color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste(model_name, "Forecast Distribution"),
         x = "Forecast Volatility", y = "Density") +
    theme_minimal()
  ggsave(paste0("figures/distribution_", model_name, ".png"), p, width = 8, height = 4, dpi = 150)
}
plot_distribution(fm_lstm, "LSTM")
plot_distribution(fm_ets, "ETS")

cat("\nDONE. All CSVs and plots saved.\n")