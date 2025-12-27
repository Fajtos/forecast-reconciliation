
data {
  int<lower=1> T;          // number of time points (origins)
  int<lower=1> H;          // forecast horizon
  int<lower=1> N;          // total number of forecasts
  int<lower=1> time_idx[N]; // target time for each forecast
  int<lower=1> horiz_idx[N]; // horizon for each forecast
  real y_forecast[N];      // AI forecasts
}
parameters {
  real<lower=-1, upper=1> phi;          // AR(1) coefficient
  real<lower=0> sigma_eta;              // state innovation SD
  vector[H] delta_raw;                  // unconstrained biases
  vector<lower=0>[H] tau;               // horizon-specific SD
  vector[T+H-1] mu;                     // latent states
}
transformed parameters {
  vector[H] delta;  // centered biases (sum to zero)
  delta = delta_raw - mean(delta_raw);
}
model {
  // Priors (as in paper)
  phi ~ normal(0.8, 0.1);
  sigma_eta ~ inv_gamma(3, 2);
  delta_raw ~ normal(0, 0.5);
  tau ~ inv_gamma(5, 2);
  
  // State evolution
  mu[1] ~ normal(0, sigma_eta/sqrt(1 - phi^2));
  for (t in 2:(T+H-1)) {
    mu[t] ~ normal(phi * mu[t-1], sigma_eta);
  }
  
  // Observation model
  for (n in 1:N) {
    y_forecast[n] ~ normal(mu[time_idx[n]] + delta[horiz_idx[n]], tau[horiz_idx[n]]);
  }
}
generated quantities {
  // For posterior predictive checks
  real y_rep[N];
  real log_lik[N];
  
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[time_idx[n]] + delta[horiz_idx[n]], tau[horiz_idx[n]]);
    log_lik[n] = normal_lpdf(y_forecast[n] | mu[time_idx[n]] + delta[horiz_idx[n]], tau[horiz_idx[n]]);
  }
}

