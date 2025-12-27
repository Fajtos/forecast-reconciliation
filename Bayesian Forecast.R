print("Most recent dates in subset:")
print(tail(dates_subset, 10))

# Use the last 30 days of our subset
window_size <- 30
window_start <- length(dates_subset) - window_size + 1
window_end <- length(dates_subset)

print(paste("\nUsing window:", window_start, "to", window_end))
print(paste("Dates:", dates_subset[window_start], "to", dates_subset[window_end]))

# Get data for this window
recent_dates <- dates_subset[window_start:window_end]
recent_original <- forecast_subset[window_start:window_end, 1]
recent_reconciled <- recon_final$mu_mean[window_start:window_end]
recent_actual <- volatility[window_start:window_end]

# Create data frame
recent_data <- data.frame(
  Date = as.Date(recent_dates),
  Original = recent_original,
  Reconciled = recent_reconciled,
  Actual = recent_actual
)

print("\nRecent data (first 5 rows):")
print(head(recent_data))

print("\nRecent data (last 5 rows):")
print(tail(recent_data))

# Create the plot with recent dates
p3_recent <- ggplot(recent_data, aes(x = Date)) +
  geom_line(aes(y = Original, color = "Original Forecast"), linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = Reconciled, color = "Bayesian Reconciled"), linewidth = 1.2) +
  geom_line(aes(y = Actual, color = "Actual Volatility"), linewidth = 0.8, linetype = "dashed") +
  labs(
    title = "Original vs Reconciled Volatility Forecasts",
    subtitle = paste("Most Recent", window_size, "Trading Days"),
    y = "Volatility",
    color = "Series"
  ) +
  scale_color_manual(
    values = c(
      "Original Forecast" = "red",
      "Bayesian Reconciled" = "blue",
      "Actual Volatility" = "black"
    ),
    labels = c(
      "Original Forecast" = "Raw AR(1) Forecast",
      "Bayesian Reconciled" = "Bayesian Reconciled",
      "Actual Volatility" = "Actual Volatility"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  ) +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d, %Y")

# Save the plot
ggsave("plots/timeseries_recent_final.png", p3_recent, width = 12, height = 6)
print("Saved: plots/timeseries_recent_final.png")

# Also update the example target plot to use recent date
# Find a recent target with high incoherence
print("\nFinding recent target with high incoherence...")

# Get targets from recent window
recent_targets <- window_start:window_end
incoherence_recent <- incoherence_detailed$details[recent_targets, ]

# Find target with highest variance (most incoherence)
high_var_target <- incoherence_recent$target[which.max(incoherence_recent$variance)]
print(paste("Recent target with highest incoherence:", high_var_target))
print(paste("Date:", target_dates[high_var_target]))
print(paste("Variance:", round(incoherence_recent$variance[which.max(incoherence_recent$variance)], 7)))

# Display the plot
print(p3_recent)