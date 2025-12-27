
data {
  int<lower=1> T;           // Number of time points
  int<lower=1> H;           // Maximum forecast horizon
  int<lower=1> N_obs;       // Number of forecast observations
  int<lower=1, upper=T> time_idx[N_obs];      // Time index for each forecast
  int<lower=1, upper=H> horizon_idx[N_obs];   // Horizon index for each forecast
  real forecast_value[N_obs];                 // Forecast values
}

parameters {
  // Latent states
  vector[T+H] mu;
  
  // AR(1) parameters
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma_eta;
  
  // Horizon-specific parameters
  vector[H] delta_raw;      // Unconstrained biases
  vector<lower=0>[H] tau;   // Horizon-specific reliability
  
  // Future: could add observation noise for y if needed
  // real<lower=0> sigma_y;
}

transformed parameters {
  // Impose sum-to-zero constraint on biases
  vector[H] delta = delta_raw - mean(delta_raw);
  
  // Likelihood contribution could be computed here if needed
}

model {
  // Priors
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
  // Posterior predictive checks
  vector[N_obs] forecast_rep;
  real log_lik[N_obs];
  
  for(n in 1:N_obs) {
    int t = time_idx[n];
    int h = horizon_idx[n];
    forecast_rep[n] = normal_rng(mu[t + h] + delta[h], tau[h]);
    log_lik[n] = normal_lpdf(forecast_value[n] | mu[t + h] + delta[h], tau[h]);
  }
}
