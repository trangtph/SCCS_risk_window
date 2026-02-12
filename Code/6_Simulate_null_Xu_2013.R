###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################


################################################################################
# Script: Generate null distributions for p-value calculation (Xu_2013) ########
################################################################################

### Load functions and packages

source(here("Code", "3_Functions_3methods.R"))
source(here("Code", "0_Helper_functions.R"))

library(dplyr)
library(SCCS)
library(survival)
library(fixest)
library(here)
library(rio)
library(foreach)


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

results_raw_3meth <- import(here("Report", "Results_raw_3methods.csv"))

### Specify parameters ---------------------------------------------------------

n_sim <- 1000
calendar_interval <- 30
exposure_scan <- c("bnt_1", "bnt_2", "chad_1", "chad_2")
output_dir <- here("Report", "Null_dist")

### Implement the analysis -----------------------------------------------------


null_dist_myo <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$myo_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  output_dir = output_dir
) 

null_dist_pe <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$pe_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  output_dir = output_dir
) 

null_dist_thrc <- loop_4_exp_sim(
  seed_list = get_seeds(n_sim = n_sim),
  n_sim = n_sim,
  data = cdm_sccs$thrc_sccs,
  expo_to_scan = exposure_scan,
  calendar_interval = calendar_interval,
  output_dir = output_dir
) 
