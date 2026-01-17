#####################################################################
## This file describe analytic code for Hunsberger scan method (2015)
## Trang's version
#####################################################################

# The method of Hunsberger (2015) scan method uses a scan statistics to identify exposure window
# H0: uniform distribution of events across the period
# H1: non-uniform distribution (clustering) of events in the period
# Steps:
#  1. Format data into compatible format  - event count by follow-up days
#  2. For each scanning window, calculate the generalized likelihood ratio test statistic (GLRT): 
#     The ratio of the maximized likelihood of events under the alternative over the max likelihood under the null 
#  3. The optimal window is the one that maximize the GLRT
#  4. Calculate the P-value for the probability of seeing the temporal cluster if events follow uniform distribution
#  5. Calculate excess risk (ER): the ratio of the estimate of the number of events
#     per day in the exposure period to the estimate of the number of events per day 
#     under the null model of events being uniformly distributed across the follow up period
#  6. Bootstrap confidence interval for the excess risk

# 1. Function for data preparation ----------
# Before this step the data should be formated to the SCCS "multi" format
prepare_scan_data <- function(data, id_col, event_col, exposure_col, max_days) {
  # Calculate event time, day 0 = day of vaccination
  data$event_time <- data[[event_col]] - data[[exposure_col]]
  data <- data[order(data[[id_col]], data$event_time), ]
  
  # Filter events within [1, max_days]
  df <- data[data$event_time >= 1 & data$event_time <= max_days, ]
  
  # Keep only the first event for each case ID
  df <- df[!duplicated(df[[id_col]]), ]
  
  # Tabulate counts up to the flexible max_days
  counts <- tabulate(df$event_time, nbins = max_days)
  
  return(list(
    event_days = df$event_time,
    daily_counts = counts,
    n = sum(counts),
    m = max_days
  ))
}

# 2. Function to calculate GLR statistic on log scale ----------
# Take the log of the GLR stat to avoid working with extremely small numbers
# log is a monotone transformation so maximizing the log is equivalent to 
# maximizing the original test stat

log_glr_stat <- function(Ec, # Number of events up to day c
                         n,  # Total number of event (= number of individuals)
                         c,  # The day being scanned
                         m) {# The maximum scanning window  
  
  if (Ec == 0 || Ec == n) return(0) # the glr stat = 1 in these cases -> log = 0
  
  p_hat <- Ec / n 
  theta <- c / m
  
  if (p_hat <= theta) return(0) # because when p_hat <= theta the glr stat = 1 -> log = 0
  
  # If p_hat > theta, return the log of the glr stat:
  Ec * log(p_hat / theta) +
    (n - Ec) * log((1 - p_hat) / (1 - theta))
}

# 3. Function to scan the GLRT over all c and choose the maximum ----------

find_max_glrt_scan <- function(counts # Daily event counts up to day m
) {
  
  n <- sum(counts)
  m <- length(counts)
  cum_events <- cumsum(counts)
  
  best_log_glr_stat <- 0
  best_c <- NA
  best_Ec <- NA
  
  for (c in 1:m) {
    # Calculate the log glr stat for each c
    Ec <- cum_events[c] # Number of events from day 1 up to day c
    llr <- log_glr_stat(Ec, n, c, m)
    
    # If this window's GLRT is higher than previous windows, save it as the 'best'
    if (llr > best_log_glr_stat) {
      best_log_glr_stat <- llr
      best_c <- c
      best_Ec <- Ec
    }
  }
  
  return(list(
    max_log_glr_stat = best_log_glr_stat,
    optimal_c = best_c,
    optimal_Ec = best_Ec
  ))
}



# 4. Function to calculate the Monte Carlo p-value -------
# for testing whether events occur earlier than would be expected 
# if events follow a uniform distribution
scan_p_value <- function(counts, n_sim = 4999, seed = 2026) {
  
  obs <- find_max_glrt_scan(counts) # The observed test stat
  m <- length(counts)
  n <- sum(counts)
  
  set.seed(seed)
  # Simulate the test stat distribution under the null
  sim_stats <- replicate(n_sim, {
    # Simulate event time for n individual under the null (uniform distribution)
    sim_event_time <- sample(1:m, size = n, replace = TRUE) 
    # Calculate the event count per day
    sim_counts <- tabulate(sim_event_time, nbins = m)
    # Calculate the test stat under the null
    find_max_glrt_scan(sim_counts)$max_log_glr_stat
  })
  # Calculate p-value: number of times the test stat under the null >= the observed one
  p_val <- (sum(sim_stats >= obs$max_log_glr_stat) + 1) / (n_sim + 1)
  
  list(
    p_value = p_val,
    observed_test_stat = obs
  )
}

# 5. Function to calculate Excess Risk ------
calculate_ER <- function(Ec, c, n, m) {
  (Ec / c) / (n / m)
}

# 6. Non-parametric bootstrap for 95% CI of ER -----
bootstrap_ER <- function(event_time, m, n_boot = 5000, seed = 2026) {
  
  n <- length(event_time)
  ERs <- numeric(n_boot) # Vector to store the bootstraped ERs
  
  # Bootstrap the ER 
  for (b in 1:n_boot) {
    # Resample the observed event time for n individuals
    boot_days <- sample(event_time, size = n, replace = TRUE)
    boot_counts <- tabulate(boot_days, nbins = m)
    
    # Recalculate the test stat and the ER for each bootstrapped sample
    res <- find_max_glrt_scan(boot_counts)
    
    ERs[b] <- calculate_ER(
      Ec = res$optimal_Ec,
      c  = res$optimal_c,
      n  = n,
      m  = m
    )
  }
  # Percentile CI
  quantile(ERs, c(0.025, 0.975), na.rm = TRUE)
}

# 7. Implement the full procedure -----
run_scan_test <- function(prepared_data,
                          n_sim = 4999,
                          n_boot = 5000, 
                          seed = 2026) {
  # Prepare data
  counts <- prepared_data$daily_counts
  n <- prepared_data$n
  m <- prepared_data$m
  event_days <- prepared_data$event_days
  
  # Scan to choose the optimal window and calculate the p-value
  scan_res <- scan_p_value(counts, n_sim, seed)
  
  # Calclate ER
  ER_hat <- calculate_ER(
    scan_res$observed$optimal_Ec,
    scan_res$observed$optimal_c,
    n, m
  )
  
  # Calculate the 95% CI of ER
  ci <- bootstrap_ER(event_time = event_days, 
                     m = m, n_boot = n_boot,
                     seed = seed)
  
  return(list(
    total_event = n,
    optimal_window_days = scan_res$observed$optimal_c,
    p_value = scan_res$p_value,
    events_in_window = scan_res$observed$optimal_Ec,
    excess_risk = ER_hat,
    er_ci_95 = ci
  ))
}

# Test the procedure --------------
## itpdat ----
itp_prepared <- prepare_scan_data(
  data = itpdat, 
  id_col = "case",
  event_col = "itp", 
  exposure_col = "mmr", 
  max_days = 147
)

itp_result_huns <- run_scan_test(prepared_data = itp_prepared)
# 17 cases in total, Optimal window: Days 1-51, p=0.0002, ER = 2.543, 95% CI [1.82, 4.32]

## gbs data ----
# Create GBS dataset
gbs_event_data <- c(1,1,2,3,3,4,6,6,7,8,9,9,10,10,11,11,12,12,12,13,14,14,16,16,17,17,17,18,19,20,
                    20,23,24,25,26,27,28,28,28,29,30,32,33,33,33,33,35,36,37,37,39,40,41,42,44,
                    46,50,50,53,54,63,64,64,67,68,71,71,74,75,76,77,84,85,87,87,87,90,90,90)

# Define study parameters
m_total <- 90  # Maximum follow-up day
n_total <- length(gbs_event_data)

# Aggregate into daily counts (Day 1, Day 2, ..., Day 90)
daily_counts <- tabulate(gbs_event_data, nbins = m_total)

# Create the prepared data list for the analysis functions
prepared_data <- list(event_days = gbs_event_data, daily_counts = daily_counts, n = n_total, m = m_total)
# Run analysis
gbs_result_huns <- run_scan_test(prepared_data = prepared_data)

