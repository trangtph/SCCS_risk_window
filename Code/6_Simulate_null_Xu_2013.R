###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################


################################################################################
# Script: Running three methods Xu_2011, Xu_2013 and Campos on the data ########
################################################################################

## -----------------------------------------------------------------------------
## Part 2: Generate null distributions for p-value calculation (Xu_2013) -------
## -----------------------------------------------------------------------------

### Specify parameters ---------------------------------------------------------

n_sim <- 2000
calendar_interval <- 30
dataset <- c("myo_sccs", "pe_sccs", "thrc_sccs") 

### Implement the analysis -----------------------------------------------------

null_dist_all_data <- foreach(dat = dataset,
                              .combine = rbind) %do% 
  {
    message(paste("Analysing dataset:", dat, "at", Sys.time()))
    
    if (!dat %in% names(cdm_sccs)) {
      stop("Dataset '", dat, "' not found in cdm_sccs")
    }
    
    data <- cdm_sccs[[dat]]
    
    sim_null_dist_xu2013(
      seed_list = get_seeds(n_sim = n_sim),
      n_sim = n_sim,
      data = data,
      calendar_interval = calendar_interval,
    ) 
  }

### Export the results ---------------------------------------------------------

write.csv(null_dist_all_data, here("Report", "Null_dist_all_data.csv"))
