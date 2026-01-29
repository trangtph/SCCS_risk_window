#####################################################################
## This file describe analytic code for Hunsberger scan method (2015)
## Trang's version
## Edited by Azida to include multiple exposure
## Last update: 28-Jan-2026
#####################################################################

# The method of Hunsberger (2015) scan method uses a scan statistics to identify exposure window
# H0: uniform distribution of events across the period
# H1: non-uniform distribution (clustering) of events in the period
# Steps:
#  1. Format data into compatible format  - event count by follow-up days
#     - define to loop by number of exposure columns in the dataset
#  2. For each scanning window, calculate the generalized likelihood ratio test statistic (GLRT): 
#     The ratio of the maximized likelihood of events under the alternative over the max likelihood under the null 
#  3. The optimal window is the one that maximize the GLRT
#  4. Calculate the P-value for the probability of seeing the temporal cluster if events follow uniform distribution
#  5. Calculate excess risk (ER): the ratio of the estimate of the number of events
#     per day in the exposure period to the estimate of the number of events per day 
#     under the null model of events being uniformly distributed across the follow up period
#  6. Bootstrap confidence interval for the excess risk


# Load necessary packages
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)
pacman::p_load(
rio, 
here
)

# ------------------------------------------------------------------------------
# Part 1: Define functions to implement the method Hunsberger (2015) -----------
# ------------------------------------------------------------------------------

## 1. Function for data preparation --------------------------------------------
# This function loops over 4 exposure columns
prepare_scan_data_multi <- function(data, id_col, event_col, exposure_cols, 
                                    start_col, end_col, max_days) {
  
  # Define dose pairs for censoring
  dose_pairs <- list(
    bnt_1  = "bnt_2",
    chad_1 = "chad_2")
  
  
  # Iterate over each exposure column provided
  results <- lapply(exposure_cols, function(exp_col) {
    
    temp_data <- data
    # Calculate event time: day 0 = day of specific exposure
    temp_data$event_time <- temp_data[[event_col]] - temp_data[[exp_col]]
    
    # Default: user-specified max_days for everyone
    temp_data$max_days_eff <- max_days
    
    # Apply dose-2 censoring if this is a first dose
    if (exp_col %in% names(dose_pairs)) {
      
      dose2_col <- dose_pairs[[exp_col]]
      has_dose2 <- !is.na(temp_data[[dose2_col]]) & 
        !is.na(temp_data[[exp_col]])
      # Calculate dosing interval
      dose_gap <- temp_data[[dose2_col]] - temp_data[[exp_col]] - 1
      # Censor max_days at dose 2
      temp_data$max_days_eff[has_dose2] <- 
        pmin(max_days, dose_gap[has_dose2])
    }
    
    # Filter conditions
    # 1. Event must be between the start and end of observation period
    # 2. event_time must be within [1, max_days]
    # 3. no missing values in critical columns (event, exposure)
         
    valid_rows <- !is.na(temp_data$event_time) & 
      !is.na(temp_data[[event_col]]) &
      temp_data[[event_col]] >= temp_data[[start_col]] & 
      temp_data[[event_col]] <= temp_data[[end_col]] &
      temp_data$event_time >= 1 & 
      temp_data$event_time <= temp_data$max_days_eff
         
    df <- temp_data[valid_rows, ]
    
    # Sort by ID and event time
    df <- df[order(df[[id_col]], df$event_time), ]
    
    # Keep only the first event for each case ID for this exposure 
    df <- df[!duplicated(df[[id_col]]), ]
    
    # Tabulate counts
    counts <- tabulate(df$event_time, nbins = max_days)
    
    return(list(
      event_days = df$event_time,
      daily_counts = counts,
      n = sum(counts),
      m = max_days
    ))
  })
  
  # Name the list elements for easy access
  names(results) <- exposure_cols
  return(results)
}

## 2. Function to calculate GLR statistic on log scale -------------------------
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

## 3. Function to scan the GLRT over all c and choose the maximum --------------

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


## 4. Function to calculate the Monte Carlo p-value ----------------------------
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

## 5. Function to calculate Excess Risk ----------------------------------------
calculate_ER <- function(Ec, c, n, m) {
  (Ec / c) / (n / m)
}

## 6. Non-parametric bootstrap for 95% CI of ER --------------------------------
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

## 7. Function to run scan test ------------------------------------------------
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
  
  # Calculate ER
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

## 8. Function to implement the full procedure for multiple exposures ----------
run_full_analysis <- function(data, id_col, event_col, exposure_cols, start_col, end_col, max_days, 
                              n_sim = 4999, n_boot = 5000, seed = 2026) {
  
  # Step A: Prepare data for all exposure columns
  prepared_list <- prepare_scan_data_multi(data, id_col, event_col, exposure_cols, 
                                           start_col, end_col, max_days)
  
  # Step B: Run the scan test on each prepared dataset
  final_results <- lapply(names(prepared_list), function(name) {
    prep <- prepared_list[[name]]
    
    # Check if there are enough events to run analysis
    if (prep$n < 2) {
      return(list(exposure = name, error = "Insufficient events"))
    }
    
    # Execute run_scan_test logic
    res <- run_scan_test(prep, n_sim = n_sim, n_boot = n_boot, seed = seed)
    res$exposure <- name
    return(res)
  })
  
  names(final_results) <- exposure_cols
  return(final_results)
}

# ------------------------------------------------------------------------------
# Part 2: Implement the method on the data -------------------------------------
# ------------------------------------------------------------------------------

## Load data -----------------------------------------------------------------

files <- c("individual", "myo_sccs", "pe_sccs", "thrc_sccs", "covid_infection", "comorb_stat", "covid_vacc_wide")

cdm_sccs <- setNames(
  lapply(files, function(nm) {
    rio::import(here("Processed_data", paste0(nm, ".rds")), trust = TRUE)
  }),
  files
)
str(cdm_sccs, max.level = 1)

## Specify parameters ----------------------------------------------------------
exposures <- c('bnt_1', 'bnt_2', 'chad_1', 'chad_2')
scan_length <- 70

## Implement the analysis ------------------------------------------------------
results <- run_full_analysis(data, "id", "myo", exposures, "obs_sta", "obs_end", scan_length)



## Check data (optional)
data_prepared <- prepare_scan_data_multi(data = cdm_eligible$pe_sccs, 
                                         id_col = "id", 
                                         event_col =  "diagnosis", 
                                         exposure_cols = exposures, 
                                         start_col = "obs_sta", end_col = "obs_end", 
                                         max_days = 185)

# Convert to a summary table
summary_df <- do.call(rbind, lapply(names(data_prepared), function(name) {
  x <- data_prepared[[name]]
  data.frame(
    exposure = name,
    total_events = x$n,
    max_days = x$m
  )
}))

print(summary_df)

daily_df <- do.call(rbind, lapply(names(data_prepared), function(name) {
  x <- data_prepared[[name]]
  data.frame(
    exposure = name,
    day = 1:x$m,
    count = x$daily_counts
  )
}))
