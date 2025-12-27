library(rstan)
library(dplyr)
library(ggplot2)

source("simulation_setup.R")

# Load the already compiled model
cat("Using pre-compiled Stan model...\n")

# Function to run one replication with better settings
run_improved_replication <- function(replication_id, T = 100, H = 5) {
  cat(sprintf("\n--- Replication %d ---\n", replication_id))
  
  # 1. Generate data
  sim_data <- simulate_forecast_matrix(T = T, H = H)
  
  # 2. Prepare Stan data
  forecast_matrix <- sim_data$forecasts
  T_actual <- sim_data$T
  H_actual <- sim_data$H
  
  forecast_values <- c()
  time_indices <- c()
  horizon_indices <- c()
  
  for(t in 1:T_actual) {
    for(h in 1:H_actual) {
      if(!is.na(forecast_matrix[t, h]) && (t + h) <= (T_actual + H_actual)) {
        forecast_values <- c(forecast_values, forecast_matrix[t, h])
        time_indices <- c(time_indices, t)
        horizon_indices <- c(horizon_indices, h)
      }
    }
  }
  
  stan_data <- list(
    T = T_actual,
    H = H_actual,
    N_obs = length(forecast_values),
    time_idx = time_indices,
    horizon_idx = horizon_indices,
    forecast_value = forecast_values
  )
  
  # 3. Fit model with better settings
  cat("Fitting model...\n")
  fit <- sampling(
    compiled_model,
    data = stan_data,
    chains = 2,           # 2 chains for better diagnostics
    iter = 1200,          # More iterations
    warmup = 600,         # More warmup
    refresh = 0,          # No progress updates (faster)
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 15
    )
  )
  
  # 4. Extract results
  params <- extract(fit)
  
  # Calculate key metrics
  mu_est <- apply(params$mu, 2, mean)
  delta_est <- apply(params$delta, 2, mean)
  tau_est <- apply(params$tau, 2, mean)
  phi_est <- mean(params$phi)
  sigma_eta_est <- mean(params$sigma_eta)
  
  # Create reconciled forecasts
  reconciled_matrix <- matrix(NA, nrow = T_actual, ncol = H_actual)
  for(t in 1:T_actual) {
    for(h in 1:H_actual) {
      if(!is.na(forecast_matrix[t, h]) && (t + h) <= length(mu_est)) {
        reconciled_matrix[t, h] <- mu_est[t + h] + delta_est[h]
      }
    }
  }
  
  # Calculate incoherence
  raw_incoherence <- calculate_incoherence(forecast_matrix)
  reconciled_incoherence <- calculate_incoherence(reconciled_matrix)
  
  # Check convergence
  rhat_values <- summary(fit)$summary[, "Rhat"]
  converged <- all(rhat_values < 1.05, na.rm = TRUE)
  
  # Calculate parameter recovery errors
  delta_rmse <- sqrt(mean((delta_est - sim_data$horizon_params$true_bias)^2))
  tau_rmse <- sqrt(mean((tau_est - sim_data$horizon_params$true_tau)^2))
  phi_error <- phi_est - 0.8  # True phi
  sigma_eta_error <- sigma_eta_est - 0.3  # True sigma_eta
  
  return(list(
    replication_id = replication_id,
    raw_incoherence = raw_incoherence,
    reconciled_incoherence = reconciled_incoherence,
    incoherence_reduction = (raw_incoherence - reconciled_incoherence) / raw_incoherence * 100,
    delta_rmse = delta_rmse,
    tau_rmse = tau_rmse,
    phi_error = phi_error,
    sigma_eta_error = sigma_eta_error,
    converged = converged,
    mu_correlation = cor(mu_est[1:T_actual], sim_data$mu, use = "complete.obs"),
    true_values = list(
      bias = sim_data$horizon_params$true_bias,
      tau = sim_data$horizon_params$true_tau
    ),
    estimates = list(
      delta = delta_est,
      tau = tau_est,
      phi = phi_est,
      sigma_eta = sigma_eta_est
    )
  ))
}

# Run a small study (5 replications for testing)
cat("=== Running Improved Simulation Study ===\n")
cat("5 replications with T=100, H=5\n")

n_replications <- 5
results <- list()

for(i in 1:n_replications) {
  result <- tryCatch({
    run_improved_replication(i, T = 100, H = 5)
  }, error = function(e) {
    cat(sprintf("Error in replication %d: %s\n", i, e$message))
    return(NULL)
  })
  
  if(!is.null(result)) {
    results[[i]] <- result
    cat(sprintf("Completed replication %d: Reduction = %.1f%%, Converged = %s\n",
                i, result$incoherence_reduction, result$converged))
  }
}

# Remove any failed replications
results <- results[!sapply(results, is.null)]

if(length(results) > 0) {
  # Create summary data frame
  summary_df <- data.frame(
    Replication = sapply(results, function(x) x$replication_id),
    Raw_Incoherence = sapply(results, function(x) x$raw_incoherence),
    Reconciled_Incoherence = sapply(results, function(x) x$reconciled_incoherence),
    Reduction_Pct = sapply(results, function(x) x$incoherence_reduction),
    Delta_RMSE = sapply(results, function(x) x$delta_rmse),
    Tau_RMSE = sapply(results, function(x) x$tau_rmse),
    Mu_Correlation = sapply(results, function(x) x$mu_correlation),
    Converged = sapply(results, function(x) x$converged)
  )
  
  cat("\n=== Summary Statistics ===\n")
  print(summary_df)
  
  # Calculate overall statistics
  overall_stats <- data.frame(
    Metric = c("Mean Raw Incoherence", "Mean Reconciled Incoherence", 
               "Mean Reduction (%)", "Mean δ_h RMSE", "Mean τ_h RMSE",
               "Mean μ Correlation", "Convergence Rate (%)"),
    Value = c(mean(summary_df$Raw_Incoherence),
              mean(summary_df$Reconciled_Incoherence),
              mean(summary_df$Reduction_Pct),
              mean(summary_df$Delta_RMSE),
              mean(summary_df$Tau_RMSE),
              mean(summary_df$Mu_Correlation, na.rm = TRUE),
              mean(summary_df$Converged) * 100)
  )
  
  cat("\n=== Overall Performance ===\n")
  print(overall_stats)
  
  # Create LaTeX table for paper
  cat("\n=== LaTeX Table for Paper ===\n")
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\begin{tabular}{@{}lcccc@{}}\n")
  cat("\\toprule\n")
  cat("\\textbf{Replication} & \\textbf{Raw $\\mathcal{I}$} & \\textbf{Reconciled $\\mathcal{I}$} & \\textbf{Reduction} & \\textbf{$\\mu$ Correlation} \\\\\n")
  cat("\\midrule\n")
  
  for(i in 1:nrow(summary_df)) {
    cat(sprintf("%d & %.4f & %.4f & %.1f\\%% & %.3f \\\\\n",
                summary_df$Replication[i],
                summary_df$Raw_Incoherence[i],
                summary_df$Reconciled_Incoherence[i],
                summary_df$Reduction_Pct[i],
                summary_df$Mu_Correlation[i]))
  }
  
  cat("\\midrule\n")
  cat(sprintf("\\textbf{Mean} & %.4f & %.4f & %.1f\\%% & %.3f \\\\\n",
              mean(summary_df$Raw_Incoherence),
              mean(summary_df$Reconciled_Incoherence),
              mean(summary_df$Reduction_Pct),
              mean(summary_df$Mu_Correlation, na.rm = TRUE)))
  cat(sprintf("\\textbf{SD} & %.4f & %.4f & %.1f\\%% & %.3f \\\\\n",
              sd(summary_df$Raw_Incoherence),
              sd(summary_df$Reconciled_Incoherence),
              sd(summary_df$Reduction_Pct),
              sd(summary_df$Mu_Correlation, na.rm = TRUE)))
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\caption{Incoherence reduction across simulation replications. Bayesian reconciliation consistently reduces temporal incoherence while accurately recovering the latent trajectory (high $\\mu$ correlation).}\n")
  cat("\\label{tab:simulation_results}\n")
  cat("\\end{table}\n")
  
  # Create visualization
  create_performance_plot <- function(summary_df) {
    # Melt data for plotting
    plot_data <- summary_df %>%
      select(Replication, Raw_Incoherence, Reconciled_Incoherence) %>%
      pivot_longer(cols = c(Raw_Incoherence, Reconciled_Incoherence),
                   names_to = "Type",
                   values_to = "Incoherence")
    
    p <- ggplot(plot_data, aes(x = factor(Replication), y = Incoherence, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("Raw_Incoherence" = "red", 
                                   "Reconciled_Incoherence" = "blue"),
                        labels = c("Raw Forecasts", "Reconciled Forecasts")) +
      labs(title = "Incoherence Reduction Across Replications",
           x = "Replication",
           y = "Incoherence (ℐ)",
           fill = "") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    
    return(p)
  }
  
  performance_plot <- create_performance_plot(summary_df)
  print(performance_plot)
  ggsave("incoherence_reduction_plot.png", performance_plot, 
         width = 10, height = 6, dpi = 300)
  
  # Create parameter recovery plot across replications
  create_parameter_recovery_summary <- function(results) {
    # Combine all parameter estimates
    all_bias <- list()
    all_tau <- list()
    
    for(res in results) {
      all_bias[[length(all_bias) + 1]] <- data.frame(
        replication = res$replication_id,
        horizon = 1:length(res$true_values$bias),
        true = res$true_values$bias,
        estimated = res$estimates$delta,
        parameter = "δ_h"
      )
      
      all_tau[[length(all_tau) + 1]] <- data.frame(
        replication = res$replication_id,
        horizon = 1:length(res$true_values$tau),
        true = res$true_values$tau,
        estimated = res$estimates$tau,
        parameter = "τ_h"
      )
    }
    
    bias_df <- do.call(rbind, all_bias)
    tau_df <- do.call(rbind, all_tau)
    combined_df <- rbind(bias_df, tau_df)
    
    # Calculate recovery statistics by horizon
    recovery_stats <- combined_df %>%
      group_by(horizon, parameter) %>%
      summarise(
        true_mean = mean(true),
        est_mean = mean(estimated),
        bias = mean(estimated - true),
        rmse = sqrt(mean((estimated - true)^2)),
        .groups = "drop"
      )
    
    # Plot
    p <- ggplot(recovery_stats, aes(x = horizon)) +
      geom_point(aes(y = true_mean, color = "True"), size = 3, shape = 17) +
      geom_point(aes(y = est_mean, color = "Estimated"), size = 3) +
      geom_errorbar(aes(ymin = est_mean - rmse, ymax = est_mean + rmse),
                    width = 0.2, color = "red", alpha = 0.5) +
      facet_wrap(~ parameter, scales = "free_y") +
      scale_color_manual(values = c("True" = "blue", "Estimated" = "red")) +
      labs(title = "Parameter Recovery Across Replications",
           x = "Horizon (h)",
           y = "Parameter Value",
           color = "",
           caption = "Error bars show ± RMSE") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    return(p)
  }
  
  param_plot <- create_parameter_recovery_summary(results)
  print(param_plot)
  ggsave("parameter_recovery_summary.png", param_plot, 
         width = 10, height = 6, dpi = 300)
  
  # Save all results
  saveRDS(list(
    results = results,
    summary = summary_df,
    overall_stats = overall_stats
  ), file = "improved_simulation_results.rds")
  
  cat("\n=== Simulation Complete ===\n")
  cat(sprintf("Successfully completed %d replications\n", length(results)))
  cat(sprintf("Mean incoherence reduction: %.1f%%\n", mean(summary_df$Reduction_Pct)))
  cat(sprintf("Mean μ correlation: %.3f\n", mean(summary_df$Mu_Correlation, na.rm = TRUE)))
  cat(sprintf("Convergence rate: %.0f%%\n", mean(summary_df$Converged) * 100))
  cat("\nFiles saved:\n")
  cat("1. improved_simulation_results.rds\n")
  cat("2. incoherence_reduction_plot.png\n")
  cat("3. parameter_recovery_summary.png\n")
  
} else {
  cat("\nNo replications completed successfully.\n")
}