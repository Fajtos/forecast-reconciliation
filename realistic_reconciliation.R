realistic_kalman_reconcile <- function(forecast_matrix, phi = 0.95, Q = 0.1) {
  T <- nrow(forecast_matrix)
  H <- ncol(forecast_matrix)
  
  # Initialize
  mu_est <- numeric(T + H)
  P <- numeric(T + H) + 1  # Start with some uncertainty
  mu_est[1] <- mean(forecast_matrix[1:min(10, T), 1], na.rm = TRUE)
  
  # Estimate measurement noise per horizon
  R <- numeric(H)
  for(h in 1:H) {
    R[h] <- var(forecast_matrix[, h], na.rm = TRUE) * 0.3  # Less aggressive
  }
  
  # Track updates per time point
  update_counts <- numeric(T + H)
  
  # First pass: collect all forecasts for each time
  for(t in 1:T) {
    for(h in 1:min(H, T + H - t)) {
      idx <- t + h
      if(!is.na(forecast_matrix[t, h])) {
        # Simple averaging for initialization
        if(update_counts[idx] == 0) {
          mu_est[idx] <- forecast_matrix[t, h]
        } else {
          mu_est[idx] <- (mu_est[idx] * update_counts[idx] + forecast_matrix[t, h]) / 
            (update_counts[idx] + 1)
        }
        update_counts[idx] <- update_counts[idx] + 1
      }
    }
  }
  
  # Second pass: smoother with uncertainty
  mu_smooth <- mu_est
  P_smooth <- P
  
  # Forward pass
  for(idx in 2:(T+H)) {
    # Predict
    mu_pred <- phi * mu_smooth[idx-1]
    P_pred <- phi^2 * P_smooth[idx-1] + Q
    
    # If we have direct forecasts for this time, update
    if(update_counts[idx] > 0) {
      # Use average of forecasts as measurement
      K <- P_pred / (P_pred + R[min(idx-1, H)]/update_counts[idx])
      measurement <- mu_est[idx]  # Our initial average
      mu_smooth[idx] <- mu_pred + K * (measurement - mu_pred)
      P_smooth[idx] <- (1 - K) * P_pred
    } else {
      mu_smooth[idx] <- mu_pred
      P_smooth[idx] <- P_pred
    }
  }
  
  # Backward pass (smoother)
  mu_final <- mu_smooth
  for(idx in (T+H-1):1) {
    C <- phi * P_smooth[idx] / (phi^2 * P_smooth[idx] + Q)
    mu_final[idx] <- mu_smooth[idx] + C * (mu_final[idx+1] - phi * mu_smooth[idx])
  }
  
  # Create reconciled forecasts with uncertainty
  reconciled <- matrix(NA, T, H)
  reconciled_sd <- matrix(NA, T, H)
  
  for(t in 1:T) {
    for(h in 1:min(H, T + H - t)) {
      idx <- t + h
      reconciled[t, h] <- mu_final[idx]
      reconciled_sd[t, h] <- sqrt(P_smooth[idx])
    }
  }
  
  return(list(
    mu_recon = mu_final,
    reconciled_forecasts = reconciled,
    reconciled_sd = reconciled_sd,
    uncertainty = sqrt(P_smooth)
  ))
}

# Test with new method
cat("=== TESTING REALISTIC RECONCILIATION ===\n")
recon1_real <- realistic_kalman_reconcile(scenario1$forecasts)
recon2_real <- realistic_kalman_reconcile(scenario2$forecasts)
recon3_real <- realistic_kalman_reconcile(scenario3$forecasts)

# Fixed coverage function that uses mean ± 1.645*sd for 90% interval
calculate_coverage_proper <- function(y_true, reconciled_list, level = 0.9) {
  reconciled <- reconciled_list$reconciled_forecasts
  reconciled_sd <- reconciled_list$reconciled_sd
  
  T <- length(y_true)
  H <- ncol(reconciled)
  z <- qnorm(1 - (1-level)/2)  # 1.645 for 90%
  
  coverage_count <- 0
  total_count <- 0
  
  for(t in (H+1):T) {
    # Use the most recent forecast for time t
    s <- t - 1  # Forecast from previous time
    h <- 1
    
    if(!is.na(reconciled[s, h]) && !is.na(reconciled_sd[s, h])) {
      lower <- reconciled[s, h] - z * reconciled_sd[s, h]
      upper <- reconciled[s, h] + z * reconciled_sd[s, h]
      
      if(y_true[t] >= lower && y_true[t] <= upper) {
        coverage_count <- coverage_count + 1
      }
      total_count <- total_count + 1
    }
  }
  
  if(total_count == 0) return(NA)
  return(coverage_count / total_count)
}

# Calculate all metrics with new method
cat("\nSCENARIO 1 (Realistic):\n")
cat("Reconciled incoherence:", round(calculate_incoherence(recon1_real$reconciled_forecasts), 4), "\n")
cat("Coverage (reconciled):", round(calculate_coverage_proper(scenario1$y, recon1_real), 3), "\n")
cat("Correlation:", round(calculate_correlation(scenario1$mu, recon1_real$mu_recon), 3), "\n")

cat("\nSCENARIO 2 (Realistic):\n")
cat("Reconciled incoherence:", round(calculate_incoherence(recon2_real$reconciled_forecasts), 4), "\n")
cat("Coverage (reconciled):", round(calculate_coverage_proper(scenario2$y, recon2_real), 3), "\n")
cat("Correlation:", round(calculate_correlation(scenario2$mu, recon2_real$mu_recon), 3), "\n")

cat("\nSCENARIO 3 (Realistic):\n")
cat("Reconciled incoherence:", round(calculate_incoherence(recon3_real$reconciled_forecasts), 4), "\n")
cat("Coverage (reconciled):", round(calculate_coverage_proper(scenario3$y, recon3_real), 3), "\n")
cat("Correlation:", round(calculate_correlation(scenario3$mu, recon3_real$mu_recon), 3), "\n")