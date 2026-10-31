data {
  int<lower=1> T;
  int<lower=1> H;
  int<lower=1> N;
  array[N] int<lower=1, upper=T> origin_idx;
  array[N] int<lower=1, upper=H> horizon_idx;
  vector[N] forecast;

  real prior_phi_mean;
  real<lower=0> prior_phi_sd;
  real<lower=0> prior_sigma_eta_shape;
  real<lower=0> prior_sigma_eta_rate;
  real<lower=0> prior_delta_sd;
  real<lower=0> prior_tau_shape;
  real<lower=0> prior_tau_rate;

  real prior_c_mean;
  real<lower=0> prior_c_sd;

  int<lower=0> T_y;
  vector[T_y] y_obs;
  array[T_y] int<lower=1, upper=T + H> y_idx;
  real<lower=0> prior_sigma_y_shape;
  real<lower=0> prior_sigma_y_rate;
}

parameters {
  real c;
  vector[T + H] z;
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma_eta_sq;
  vector[H] delta_raw;
  vector<lower=0>[H] tau_sq;
  real<lower=0> sigma_y_sq;
}

transformed parameters {
  vector[H] delta = delta_raw - mean(delta_raw);
  real<lower=0> sigma_eta = sqrt(sigma_eta_sq);
  vector<lower=0>[H] tau = sqrt(tau_sq);
  real<lower=0> sigma_y = sqrt(sigma_y_sq);

  vector[T + H] mu;
  mu[1] = c + sigma_eta / sqrt(1 - square(phi)) * z[1];
  for (t in 2:(T + H)) {
    mu[t] = c + phi * (mu[t - 1] - c) + sigma_eta * z[t];
  }
}

model {
  c ~ normal(prior_c_mean, prior_c_sd);
  z ~ normal(0, 1);
  phi ~ normal(prior_phi_mean, prior_phi_sd);
  sigma_eta_sq ~ inv_gamma(prior_sigma_eta_shape, prior_sigma_eta_rate);
  delta_raw ~ normal(0, prior_delta_sd);
  tau_sq ~ inv_gamma(prior_tau_shape, prior_tau_rate);
  sigma_y_sq ~ inv_gamma(prior_sigma_y_shape, prior_sigma_y_rate);

  for (n in 1:N) {
    forecast[n] ~ normal(
      mu[origin_idx[n] + horizon_idx[n]] + delta[horizon_idx[n]],
      tau[horizon_idx[n]]
    );
  }

  if (T_y > 0) {
    y_obs ~ normal(mu[y_idx], sigma_y);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[T_y] y_rep;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(
      forecast[n] |
      mu[origin_idx[n] + horizon_idx[n]] + delta[horizon_idx[n]],
      tau[horizon_idx[n]]
    );
  }
  for (i in 1:T_y) {
    y_rep[i] = normal_rng(mu[y_idx[i]], sigma_y);
  }
}