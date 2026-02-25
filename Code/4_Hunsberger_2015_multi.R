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
library(here)

################################################################################
# Part 1: Define functions to implement the method Hunsberger (2015) -----------
################################################################################

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
    
    # Restrict to individuals who received this exposure within the observation period
    exposed <- !is.na(temp_data[[exp_col]]) &
      temp_data[[exp_col]] >= temp_data[[start_col]] &
      temp_data[[exp_col]] <= temp_data[[end_col]]
    
    temp_data <- temp_data[exposed, ]
    
    # Stop if filtered data are empty
    if (nrow(temp_data) == 0) {
      return(list(
        event_days   = numeric(0),
        daily_counts = integer(max_days),
        risk_set     = integer(max_days),
        n            = 0,
        m            = max_days
      ))
    }
    
    # Calculate event time: day 0 = day of specific exposure
    temp_data$event_time <- temp_data[[event_col]] - temp_data[[exp_col]]
    
    # Censor max_days at end_col, or for dose 1: min(end_col, dose-2)
    
    censor_date <- temp_data[[end_col]]
    
    if (exp_col %in% names(dose_pairs)) {
      dose2_col <- dose_pairs[[exp_col]]
      has_dose2 <- !is.na(temp_data[[dose2_col]])
      
      censor_date[has_dose2] <-
        pmin(censor_date[has_dose2], temp_data[[dose2_col]][has_dose2] - 1)
    }
    
    # Individual-specific max follow-up (days since exposure)
    temp_data$max_days_eff <- pmin(
      max_days,
      censor_date - temp_data[[exp_col]]
    )
    
    # In case of non-positive max_days_eff, set it to 0
    temp_data$max_days_eff[temp_data$max_days_eff < 0] <- 0
    
    # Construct the risk set: number exposed individuals who are not censored at day t
    risk_set <- tabulate(
      unlist(lapply(temp_data$max_days_eff, seq_len)),
      nbins = max_days
    ) 
    
    # Filter events: 
    # 1. Event must be between the start and end of observation period
    # 2. event_time must be within [1, max_days]
    # 3. no missing values in event, exposure
         
    valid_rows <- !is.na(temp_data[[event_col]]) &
      temp_data[[event_col]] >= temp_data[[start_col]] & 
      temp_data[[event_col]] <= temp_data[[end_col]] &
      temp_data$event_time >= 1 & 
      temp_data$event_time <= temp_data$max_days_eff
         
    df <- temp_data[valid_rows, ]
    
    # Keep only the first event for each case ID for this exposure (we do not need this because we already select first event in data cleaning)
    #df <- df[order(df[[id_col]], df$event_time), ]
    #df <- df[!duplicated(df[[id_col]]), ]
    
    # Tabulate counts
    counts <- tabulate(df$event_time, nbins = max_days)
    
    return(list(
      event_days = df$event_time,
      daily_counts = counts,
      risk_set = risk_set,
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
scan_p_value <- function(counts, risk_set, n_sim = 4999, seed = 2026) {
  
  obs <- find_max_glrt_scan(counts) # The observed test stat
  m <- length(counts)
  n <- sum(counts)
  
  # Calculate the probability of event being allocated to day t:
  # Under the null, event hazard is constant over observed person-time, conditional on exposure
  # In other words, probability of event occurring on day t is proportional to how much person-time exists on day t.
  prob <- risk_set / sum(risk_set)
  
  # Simulate the test stat distribution under the null
  set.seed(seed)
  sim_stats <- replicate(n_sim, {
    # simulate counts per day using a multinomial distribution with weights equal to daily risk set sizes
    sim_counts <- as.vector(rmultinom(1, size = n, prob = prob))   
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
bootstrap_ER <- function(event_time, m, n_boot = 5000, seed = 2026, min_events = 5) {
  
  n <- length(event_time)
  
  # Only calculate 95% CI if there are at least 5 events in the risk window 
  if (n < min_events) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  ERs <- numeric(n_boot) # Vector to store the bootstraped ERs
  
  set.seed(seed)
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
  risk_set <- prepared_data$risk_set
  n <- prepared_data$n
  m <- prepared_data$m
  event_days <- prepared_data$event_days
  
  # Scan to choose the optimal window and calculate the p-value
  scan_res <- scan_p_value(counts = counts, risk_set = risk_set, 
                           n_sim = n_sim, seed = seed)
  
  # Calculate ER
  ER_hat <- calculate_ER(
    Ec = scan_res$observed$optimal_Ec,
    c = scan_res$observed$optimal_c,
    n = n, m = m
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
    res <- run_scan_test(prepared_data = prep, n_sim = n_sim, n_boot = n_boot, seed = seed)
    res$exposure <- name
    return(res)
  })
  
  names(final_results) <- exposure_cols
  return(final_results)
}

## 9. Functions to make result table -------------------------------------------

# Only extract the field if it exists
get_or_na <- function(x, name) {
  if (!is.null(x[[name]])) x[[name]] else NA
}

results_to_table <- function(result_list) {
  
  rows <- list()
  
  for (outcome in names(result_list)) {
    for (exposure in names(result_list[[outcome]])) {
      
      res <- result_list[[outcome]][[exposure]]
      
      er_ci_low <- NA
      er_ci_up  <- NA
      
      if (!is.null(res$er_ci_95)) {
        er_ci_low <- unname(res$er_ci_95["2.5%"])
        er_ci_up  <- unname(res$er_ci_95["97.5%"])
      }
      #append the data frame as a new element at the end of the list
      rows[[length(rows) + 1]] <- data.frame(
        outcome              = outcome,
        exposure             = exposure,
        error                = get_or_na(res, "error"),
        total_event          = get_or_na(res, "total_event"),
        optimal_window_days  = get_or_na(res, "optimal_window_days"),
        p_value              = get_or_na(res, "p_value"),
        events_in_window     = get_or_na(res, "events_in_window"),
        excess_risk          = get_or_na(res, "excess_risk"),
        er_ci_low            = er_ci_low,
        er_ci_up             = er_ci_up,
        stringsAsFactors     = FALSE
      )
    }
  }
  # combining the list of dataframes
  do.call(rbind, rows)
}


################################################################################
# Part 2: Implement the method on the data -------------------------------------
################################################################################

## Load data -----------------------------------------------------------------

files <- c("individual", "myo_sccs", "pe_sccs", "thrc_sccs", "covid_infection", "comorb_stat", "covid_vacc_wide")

cdm_sccs <- setNames(
  lapply(files, function(nm) {
    readRDS(here("Processed_data", paste0(nm, ".rds")))
  }),
  files
)
str(cdm_sccs, max.level = 1)

## Specify parameters ----------------------------------------------------------
exposures <- c('bnt_1', 'bnt_2', 'chad_1', 'chad_2')
scan_length <- 70

## Implement the analysis ------------------------------------------------------
results_huns <- list(myo = NULL, pe = NULL, thrc = NULL)

results_huns[c("myo", "pe", "thrc")] <- lapply(cdm_sccs[c("myo_sccs", "pe_sccs", "thrc_sccs")], function(df) { 
  run_full_analysis(data = df, 
                    id_col = "id", event_col = "diagnosis", 
                    exposure_cols = exposures, 
                    start_col = "obs_sta", end_col = "obs_end", max_days = scan_length)
  })

results_huns_tab <- results_to_table(results_huns)

## Export the results

write.csv(results_huns_tab, here("Report", "Results_Hunsberger.csv"))


