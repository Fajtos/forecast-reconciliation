# simulation setup code

simulate_forecast_matrix <- function(T = 100, H = 5,
                                     phi = 0.8, sigma_eta = 0.3, sigma_y = 0.15,
                                     seed = NULL) {
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
  
  list(T = T, H = H, mu = mu[1:T],
       forecasts = forecast_matrix, y = y,
       true_bias = delta, true_tau = tau)
}

simulate_forecast_matrix_tar <- function(T = 100, H = 5,
                                         phi = 0.8, sigma_eta = 0.3, sigma_y = 0.15,
                                         tar_strength = 2.0,
                                         seed = NULL) {
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
  
  list(T = T, H = H, mu = mu[1:T],
       forecasts = forecast_matrix, y = y,
       true_bias = delta, true_tau = tau)
}

calculate_incoherence <- function(forecast_matrix) {
  T <- nrow(forecast_matrix); H <- ncol(forecast_matrix)
  total_sse <- 0; count <- 0
  
  # forecast from origin `i` targets `s` if `i + h = s`
  for (s in seq_len(T)) {
    vals <- numeric(0)
    for (i in max(1, s - H):(s - 1)) {
      h <- s - i
      if (h >= 1 && h <= H && !is.na(forecast_matrix[i, h])) {
        vals <- c(vals, forecast_matrix[i, h])
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