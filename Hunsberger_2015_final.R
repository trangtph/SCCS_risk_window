#####################################################################
## This file describe analytic code for Hunsberger scan method (2015)
## Last update: 2025-12
#####################################################################

# The method of Hunsberger (2015) scan method uses a scan statistics to identify exposure window
# H0: uniform distribution of events across the period
# H1: non-uniform distribution (clustering) of events across the period
# Steps:
#  1. Format data into compatible format  - event count by follow-up days
#  2. Calculate likelihood of events under (i) null and (ii) alternative hypothesis
#  3. Scan through every possible window, calculate and find max LR to get optimal window (c')
#  4. P-value for the probability of seeing the temporal cluster if events follow uniform distribution
#  5. Calculate excess risk - ratio of event rate inside the optimal window compared the event rate outside it
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


##--- Likelihood Function ---##

# Calculates likelihood for the Alternative Hypothesis (there IS a cluster)
calculate_likelihood_alt <- function(Ec,      # Nr of ppl with events within window 1-c
                                     n, 
                                     c_prime, 
                                     m) {     # max window to scan
  if (n == 0 || m == 0) return(NA) # safety check to prevent errors in mathematical calculations
  p_hat <- Ec / n        # Proportion of events that fell inside the window
  pi_hat <- c_prime / m  # Proportion of time the window represents
  # If the rate inside is higher than expected, calculate the binomial probability
  if (p_hat > pi_hat) return((p_hat^Ec) * ((1 - p_hat)^(n - Ec))) 
  # Otherwise, treat it as the null (no cluster found)
  else return((pi_hat^Ec) * ((1 - pi_hat)^(n - Ec)))
}

# Calculates likelihood for the Null Hypothesis (events are UNIFORM)
# Under the null, event are uniform on [1, m].
# Denote Yi the indicator if the event of individual i happened by day c. 
# Then under the null Yi ~ Bernoulli(theta = c/m).
# Ec = sum of n Yi -> Ec ~ Binomial(n, theta)
calculate_likelihood_null <- function(Ec, n, c_prime, m) {
  if (m == 0) return(NA) # safety check 
  pi_hat <- c_prime / m  # Probability based strictly on duration
  return((pi_hat^Ec) * ((1 - pi_hat)^(n - Ec)))
}

##--- Finding the Cluster (Scanning) ---##

# This function scans through every possible window (1 to 'c' days) to find the 'best' cluster
find_max_glrt_scan <- function(counts) {
  n <- sum(counts); m <- length(counts)
  if (n == 0 || m == 0) return(list(stat_lambda = NA, optimal_duration_c = NA, optimal_Ec = NA)) # safety check
  max_lambda <- -Inf; best_c <- NA; optimal_Ec <- NA
  cum_events <- cumsum(counts) # Speed optimization: pre-calculate running total of events
  
  for (c in 1:m) {
    E_c <- cum_events[c] # Events found from Day 1 to Day 'c'
    L_alt <- calculate_likelihood_alt(E_c, n, c, m)
    L_null <- calculate_likelihood_null(E_c, n, c, m)
    
    # Calculate the Likelihood Ratio (Lambda)
    lambda_curr <- L_alt / L_null
    
    # If this window's ratio is higher than previous windows, save it as the 'best'
    if (lambda_curr > max_lambda) { 
      max_lambda <- lambda_curr; best_c <- c; optimal_Ec <- E_c 
    }
  }
  return(list(stat_lambda = max_lambda, optimal_duration_c = best_c, optimal_Ec = optimal_Ec))
}

##--- Risk Calculation ---##
# ER = (Events in Window / Days in Window) / (Total Events / Total Days)
# Relative Risk (RR): Rate inside window divided by rate outside (or average)
calculate_excess_risk <- function(Ec, c_prime, n, m) {
  if (is.na(c_prime) || c_prime == 0 || n == 0) return(NA)
  rate_in_window <- Ec / c_prime
  overall_rate <- n / m
  return(rate_in_window / overall_rate) # Returns Relative Risk (RR)
}

# --- Simulation Logic (P-value and CI) ---
run_full_analysis <- function(prepared_data, n_sim = 999, n_boot = 499) {
  counts <- prepared_data$daily_counts
  n <- prepared_data$n
  m <- prepared_data$m
  
  # 1. Observed Result
  obs_res <- find_max_glrt_scan(counts)
  
  # Calculate Excess Risk (RR) for the observed window
  observed_ER <- calculate_excess_risk(obs_res$optimal_Ec, obs_res$optimal_duration_c, n, m)
  
  # 2. Monte Carlo P-value
  sim_stats <- replicate(n_sim, {
    sim_counts <- tabulate(sample(1:m, size = n, replace = TRUE), nbins = m)
    find_max_glrt_scan(sim_counts)$stat_lambda
  })
  p_val <- (sum(sim_stats >= obs_res$stat_lambda, na.rm = TRUE) + 1) / (n_sim + 1)
  

  # 3. Bootstrap CI for Excess Risk
  boot_ers <- replicate(n_boot, {
    sim_counts <- tabulate(sample(1:m, size = n, replace = TRUE), nbins = m)
    res <- find_max_glrt_scan(sim_counts)
    calculate_excess_risk(res$optimal_Ec, res$optimal_duration_c, n, m)
  })
  ci <- quantile(boot_ers, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return combined results
  return(list(
    p_value = p_val, 
    window_days = obs_res$optimal_duration_c, 
    events_in_window = obs_res$optimal_Ec,
    excess_risk = observed_ER,
    ci_95 = ci
  ))
}




##--- Execution ---##

### Example 1: Use 'itpdat' from SCCS package ###
library(SCCS)
data(itpdat)

# Define your parameters
case <- "case"      # Individual identifier
event <- "itp"      # Adverse event age
vaccine <- "mmr"    # Exposure age
scan_length <- 147  # Potential window length to be scanned

# Step 1: Prepare data using all flexible arguments
itp_prepared <- prepare_scan_data(
  data = itpdat, 
  id_col = "case",
  event_col = "itp", 
  exposure_col = "mmr", 
  max_days = 147
)

# Step 2: Run the scan statistic (using the previously defined run_full_analysis function)
final_results <- run_full_analysis(itp_prepared, n_sim = 999, n_boot = 499)
# 17 cases in total, Optimal window: Days 1-51, p=0.01, ER = 2.543, 95% CI [1.014, 8.647]




### Example 2: GBS Dataset from example described in the paper ###
##-- 1. Data Preparation

# Create GBS dataset
gbs_event_data <- c(1,1,2,3,3,4,6,6,7,8,9,9,10,10,11,11,12,12,12,13,14,14,16,16,17,17,17,18,19,20,
                20,23,24,25,26,27,28,28,28,29,30,32,33,33,33,33,35,36,37,37,39,40,41,42,44,
                46,50,50,53,54,63,64,64,67,68,71,71,74,75,76,77,84,85,87,87,87,90,90,90)

# 1. Define study parameters
m_total <- 90  # Maximum follow-up day
n_total <- length(gbs_event_data)

# 2. Aggregate into daily counts (Day 1, Day 2, ..., Day 90)
daily_counts <- tabulate(gbs_event_data, nbins = m_total)

# 3. Create the prepared data list for the analysis functions
prepared_data <- list(event_days = gbs_event_data, daily_counts = daily_counts, n = n_total, m = m_total)

# Run the analysis
results <- run_full_analysis(prepared_data, n_sim = 9999, n_boot = 9999)


# Format and print the final analysis summary
cat("\n===========================================")
cat("\n   TEMPORAL SCAN STATISTIC RESULTS         ")
cat("\n===========================================")
cat("\nFollow-up Duration  :", m_total, "days")
cat("\nTotal Events        :", n_total)
cat("\n-------------------------------------------")
cat("\nDetected Risk Window: Days 1 to", results$window_days)
cat("\nEvents in Window    :", results$events_in_window)
cat("\nMonte Carlo P-value :", round(results$p_value, 4))
cat("\n-------------------------------------------")
cat("\nExcess Risk (RR)    :", round(results$excess_risk, 3))
cat("\n95% Bootstrap CI    : [", round(results$ci_95[1], 3), 
    ",", round(results$ci_95[2], 3), "]")
cat("\n===========================================\n")


## Interpretation: detected optimal window = 37 days; ER 1.54; pvalue<0.05; CI 1.01, 3.42
## Similar to values reported in Hunsberger paper; p-value & CI not exactly the same due to simulation