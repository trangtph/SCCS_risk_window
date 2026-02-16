###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################


################################################################################
# Script: Running three methods Xu_2011, Xu_2013 and Campos on the data ########
################################################################################

## -----------------------------------------------------------------------------
## Part 1: Generate raw results for three methods ------------------------------
## -----------------------------------------------------------------------------


### Load functions and packages
library(dplyr)
library(SCCS)
library(survival)
library(fixest)
library(here)
library(rio)
library(foreach)
library(ggplot2)
library(readr)

source(here("Code", "3_Functions_3methods.R"))
source(here("Code", "0_Helper_functions.R"))


### Load data -------------------------------------------------------------------

files <- c("individual", "myo_sccs", "pe_sccs", "thrc_sccs", "covid_infection", "comorb_stat", "covid_vacc_wide")

cdm_sccs <- setNames(
  lapply(files, function(nm) {
    readRDS(here("Processed_data", paste0(nm, ".rds")))
  }),
  files
)
str(cdm_sccs, max.level = 1)

# make sure the id is numeric
all_ids <- unique(unlist(lapply(cdm_sccs, `[[`, "id")))

cdm_sccs <- lapply(cdm_sccs, function(df) {
  df$id_num <- as.numeric(factor(df$id, levels = all_ids))
  df
})


### Specify parameters ---------------------------------------------------------

risk_win <- seq(1, 70, by = 1)
calendar_interval <- 30
dataset <- c("myo_sccs", "pe_sccs", "thrc_sccs")

### Scanning the risk window for all exposures in three datasets ---------------

results_raw_3meth <- foreach(dat = dataset,
                             .combine = rbind) %do% 
  {
    message(paste("Analysing dataset:", dat, "at", Sys.time(), "----------------"))
    
    if (!dat %in% names(cdm_sccs)) {
      stop("Dataset '", dat, "' not found in cdm_sccs")
    }
    
    data <- cdm_sccs[[dat]]
    
    loop_4_exp(
      data = data,
      risk_win = risk_win,
      calendar_interval = 30,
      output_dir = here("Report", "Raw_results")
    )
    
  }

readr::write_csv(
  results_raw_3meth,
  here("Report", "Raw_results", "Results_raw_3methods.csv")
)

## -----------------------------------------------------------------------------
## Part 2: Obtain the optimal risk window by each method -----------------------
## -----------------------------------------------------------------------------

optimal_risk_win_3meth <- select_risk_win(data = results_raw_3meth,
                                          summary_dir = here("Report", "Summary"),
                                          plot_dir = here("Report", "Plot"),
                                          summary_file_name = "Summary_all_scens")


## -----------------------------------------------------------------------------
## Part 3: p-value calculation (Xu_2013) ---------------------------------------
## -----------------------------------------------------------------------------

### Specify parameters ---------------------------------------------------------

n_sim <- 1000
calendar_interval <- 30
exposure_scan <- c("bnt_1", "bnt_2", "chad_1", "chad_2")
output_dir <- here("Report", "Null_dist")

### Generate null distributions of the test statistic --------------------------

# Warning: it could take one to two hour to simulate data for each outcome
null_dist_myo <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$myo_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  result_table = results_raw_3meth,
  output_dir = output_dir
) 

null_dist_pe <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$pe_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  result_table = results_raw_3meth,
  output_dir = output_dir
) 

null_dist_thrc <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$thrc_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  result_table = results_raw_3meth,
  output_dir = output_dir
)

### Calculate p-value for testing the null hypothesis --------------------------
# that there does not exist an interval with elevated risk

pval_xu2013 <- p_val_cal_xu2013()


