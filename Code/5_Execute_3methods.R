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
    
    loop_risk_win(risk_win = risk_win,
                  data = data,
                  calendar_interval = calendar_interval
    )
    
  }

### Export the results ---------------------------------------------------------

write.csv(results_raw_3meth, here("Report", "Results_raw_3methods.csv"))

