library(rstan)
library(bayesplot)

test_rep <- sim_results[[1]]
T <- test_rep$T
H <- test_rep$H
forecast_matrix <- test_rep$forecast_matrix

cat("Testing Bayesian reconciliation on replication 1\n")
cat(sprintf("T = %d, H = %d\n", T, H))
cat(sprintf("Raw incoherence: %.6f\n", test_rep$raw_incoherence))

# Prepare data for Stan
prepare_stan_data <- function(forecast_matrix) {
  T <- nrow(forecast_matrix)
  H <- ncol(forecast_matrix)
  

  forecast_values <- c()
  time_indices <- c()
  horizon_indices <- c()
  
  for(t in 1:T) {
    for(h in 1:H) {
      if(!is.na(forecast_matrix[t, h]) && (t + h) <= (T + H)) {
        forecast_values <- c(forecast_values, forecast_matrix[t, h])
        time_indices <- c(time_indices, t)
        horizon_indices <- c(horizon_indices, h)
      }
    }
  }
  
  return(list(
    T = T,
    H = H,
    N_obs = length(forecast_values),
    time_idx = time_indices,
    horizon_idx = horizon_indices,
    forecast_value = forecast_values
  ))
}

stan_data <- prepare_stan_data(forecast_matrix)
cat(sprintf("Number of forecast observations: %d\n", stan_data$N_obs))

# Create Stan model
stan_model_code <- '
data {
  int<lower=1> T;           // Number of time points
  int<lower=1> H;           // Maximum forecast horizon
  int<lower=1> N_obs;       // Number of forecast observations
  int<lower=1, upper=T> time_idx[N_obs];      // Time index for each forecast
  int<lower=1, upper=H> horizon_idx[N_obs];   // Horizon index for each forecast
  real forecast_value[N_obs];                 // Forecast values
}

parameters {
  // Latent states (we only need up to T+H for all forecasts)
  vector[T+H] mu;
  
  // AR(1) parameters
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma_eta;
  
  // Horizon-specific parameters
  vector[H] delta_raw;      // Unconstrained biases
  vector<lower=0>[H] tau;   // Horizon-specific reliability
}

transformed parameters {
  // Impose sum-to-zero constraint on biases
  vector[H] delta = delta_raw - mean(delta_raw);
}

model {
  // Priors (weakly informative)
  phi ~ normal(0.8, 0.1);
  sigma_eta ~ inv_gamma(3, 2);
  delta_raw ~ normal(0, 0.5);
  tau ~ inv_gamma(5, 2);
  
  // Latent state evolution (AR(1))
  mu[1] ~ normal(0, sigma_eta/sqrt(1 - phi^2));
  for(t in 2:(T+H)) {
    mu[t] ~ normal(phi * mu[t-1], sigma_eta);
  }
  
  // Observation model for forecasts
  for(n in 1:N_obs) {
    int t = time_idx[n];
    int h = horizon_idx[n];
    forecast_value[n] ~ normal(mu[t + h] + delta[h], tau[h]);
  }
}

generated quantities {
  // For posterior predictive checks
  vector[N_obs] forecast_rep;
  real log_lik[N_obs];
  
  for(n in 1:N_obs) {
    int t = time_idx[n];
    int h = horizon_idx[n];
    forecast_rep[n] = normal_rng(mu[t + h] + delta[h], tau[h]);
    log_lik[n] = normal_lpdf(forecast_value[n] | mu[t + h] + delta[h], tau[h]);
  }
}
'

cat("Compiling Stan model...\n")
stan_model <- stan_model(model_code = stan_model_code, model_name = "forecast_reconciliation")

cat("Fitting model (this may take a few minutes)...\n")
fit <- sampling(
  stan_model,
  data = stan_data,
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("Model fitting complete!\n")

# Check convergence
cat("\n=== Convergence Diagnostics ===\n")
print(fit, pars = c("phi", "sigma_eta"))

# Extract parameters
params <- extract(fit)

# Create parameter recovery plot
plot_parameter_recovery <- function(fit, true_bias, true_tau) {
  # Extract posterior samples for delta and tau
  delta_samples <- as.matrix(fit, pars = paste0("delta[", 1:length(true_bias), "]"))
  tau_samples <- as.matrix(fit, pars = paste0("tau[", 1:length(true_tau), "]"))
  
  # Calculate posterior means and credible intervals
  delta_means <- apply(delta_samples, 2, mean)
  delta_lower <- apply(delta_samples, 2, quantile, probs = 0.05)
  delta_upper <- apply(delta_samples, 2, quantile, probs = 0.95)
  
  tau_means <- apply(tau_samples, 2, mean)
  tau_lower <- apply(tau_samples, 2, quantile, probs = 0.05)
  tau_upper <- apply(tau_samples, 2, quantile, probs = 0.95)
  
  # Create data frames for plotting
  delta_df <- data.frame(
    horizon = 1:length(true_bias),
    true = true_bias,
    estimated = delta_means,
    lower = delta_lower,
    upper = delta_upper
  )
  
  tau_df <- data.frame(
    horizon = 1:length(true_tau),
    true = true_tau,
    estimated = tau_means,
    lower = tau_lower,
    upper = tau_upper
  )
  
  # Plot delta recovery
  p1 <- ggplot(delta_df, aes(x = horizon)) +
    geom_point(aes(y = true), color = "blue", size = 3, shape = 17) +
    geom_point(aes(y = estimated), color = "red", size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "red", alpha = 0.5) +
    geom_segment(aes(x = horizon, xend = horizon, y = true, yend = estimated),
                 alpha = 0.3) +
    labs(title = "Bias Parameter (δ_h) Recovery",
         x = "Horizon (h)",
         y = "Bias",
         caption = "Blue triangles: True values\nRed points: Posterior means with 90% credible intervals") +
    theme_minimal()
  
  # Plot tau recovery
  p2 <- ggplot(tau_df, aes(x = horizon)) +
    geom_point(aes(y = true), color = "blue", size = 3, shape = 17) +
    geom_point(aes(y = estimated), color = "red", size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "red", alpha = 0.5) +
    geom_segment(aes(x = horizon, xend = horizon, y = true, yend = estimated),
                 alpha = 0.3) +
    labs(title = "Reliability Parameter (τ_h) Recovery",
         x = "Horizon (h)",
         y = "τ_h",
         caption = "Blue triangles: True values\nRed points: Posterior means with 90% credible intervals") +
    theme_minimal()
  
  return(list(delta_plot = p1, tau_plot = p2))
}

# Create parameter recovery plots
cat("\nCreating parameter recovery plots...\n")
recovery_plots <- plot_parameter_recovery(fit, 
                                          test_rep$true_bias, 
                                          test_rep$true_tau)

print(recovery_plots$delta_plot)
print(recovery_plots$tau_plot)

# Save the plots
ggsave("delta_recovery_test.png", recovery_plots$delta_plot, 
       width = 8, height = 6, dpi = 300)
ggsave("tau_recovery_test.png", recovery_plots$tau_plot, 
       width = 8, height = 6, dpi = 300)

cat("\nParameter recovery plots saved:\n")
cat("1. delta_recovery_test.png\n")
cat("2. tau_recovery_test.png\n")

# Calculate reconciled incoherence
calculate_reconciled_incoherence <- function(fit, forecast_matrix) {
  # Extract posterior mean of mu
  mu_samples <- as.matrix(fit, pars = "mu")
  mu_mean <- apply(mu_samples, 2, mean)
  
  T <- nrow(forecast_matrix)
  H <- ncol(forecast_matrix)
  
  # Create reconciled forecast matrix using posterior mean of mu
  reconciled_matrix <- matrix(NA, nrow = T, ncol = H)
  
  for(t in 1:T) {
    for(h in 1:H) {
      if(!is.na(forecast_matrix[t, h]) && (t + h) <= length(mu_mean)) {
        # For reconciliation, we use the estimated mu plus estimated delta
        # But for simple incoherence calculation, we can just use mu
        reconciled_matrix[t, h] <- mu_mean[t + h]
      }
    }
  }
  
  # Calculate incoherence of reconciled forecasts
  reconciled_incoherence <- calculate_incoherence(reconciled_matrix)
  
  return(list(
    reconciled_incoherence = reconciled_incoherence,
    reconciled_matrix = reconciled_matrix,
    mu_estimate = mu_mean
  ))
}

# Calculate reconciled incoherence
reconciliation_results <- calculate_reconciled_incoherence(fit, forecast_matrix)

cat("\n=== Reconciliation Results ===\n")
cat(sprintf("Raw incoherence: %.6f\n", test_rep$raw_incoherence))
cat(sprintf("Reconciled incoherence: %.6f\n", reconciliation_results$reconciled_incoherence))
cat(sprintf("Reduction: %.2f%%\n", 
            (test_rep$raw_incoherence - reconciliation_results$reconciled_incoherence) / 
              test_rep$raw_incoherence * 100))