###########################################
### Project: SCCS Risk Window Scanning ####
### Author: Trang - Azida              ####
###########################################

################################################################################
# Script: Defining functions for methods Xu_2011, Xu_2013 and Campos ###########
################################################################################

# ------------------------------------------------------------------------------
# Part 1: General functions to for three methods -------------------------------
# ------------------------------------------------------------------------------

## 0. Function to edit SCCS data to take into account precedence rule ----------

# Helper function to apply precedence logic for dose 2:
# If pre-exposure for dose 2 overlaps with risk period of dose 1, dose 1 takes priority
# Case 1: If gap between dose1 & dose2 <= entire risk period, pre-exposure for dose2 = NA 
# Case 2: If risk period dose 1 + pre_ex > gap between dose 1 & dose 2 > risk period dose 1, 
# the remaining days between dose1 +risk period and dose2 = pre-exposure
# Else: pre-exposure starts pre_ex days before dose 2

compute_pre_ex_dose2 <- function(dose1, dose2, risk1, pre_ex) {
  gap <- dose2 - dose1
  case_when(
    is.na(dose2) ~ NA_real_,
    gap <= risk1 + 1 ~ NA_real_,
    gap <= risk1 + pre_ex + 1 ~ dose1 + risk1 + 1,
    TRUE ~ dose2 - pre_ex
  )
}

# Function to prep the dataset:
# - Compute pre-exposure period for all exposures
# - Flexibly assign the risk window to each exposure depending on which exposure is under scanning
# - Obtain correct variable names for extracting the SCCS model results

sccs_prep <- function(data,
                      risk_scan, # The risk window we are scanning for the expo_of_interest
                      expo_of_interest, # The exposure we are scanning the risk window for,
                      pre_ex = 30, # Length of pre-exposure period
                      max_risk_win = 70, # The default risk window length for all exposures
                      calendar_adj # the calendar time groups for calendar time adjustment
                       ) {
  
  # Assign correct risk window to each exposure
  # window of exposure under scanning = risk_scan, all other exposures' windows = 70
  risk <- c(
    scan   = risk_scan,
    bnt_1  = max_risk_win,
    bnt_2  = max_risk_win,
    chad_1 = max_risk_win,
    chad_2 = max_risk_win
  )
  
  risk[expo_of_interest] <- risk["scan"]
  
  # For the null model for Xu_2013: set risk_scan = 0
  
  risk_null <- c(
    bnt_1  = max_risk_win,
    bnt_2  = max_risk_win,
    chad_1 = max_risk_win,
    chad_2 = max_risk_win
  )
  
  risk_null[expo_of_interest] <- 0
  
  
  # Compute pre-exposure period for all exposures
  data_out <- data %>%
    mutate(
      # Dose 1
      pre_bnt_1  = bnt_1  - pre_ex,
      pre_chad_1 = chad_1 - pre_ex,
      
      # Dose 2
      pre_bnt_2 = compute_pre_ex_dose2(
        dose1 = bnt_1, dose2 = bnt_2,
        risk1 = risk[["bnt_1"]], pre_ex = pre_ex
      ),
      pre_chad_2 = compute_pre_ex_dose2(
        dose1 = chad_1, dose2 = chad_2,
        risk1 = risk[["chad_1"]], pre_ex = pre_ex
      )
    )
  
  # Obtain variable names for fitting model in the later steps
  vac_interest <- sub("_.*$", "", expo_of_interest)
  dose_interest <- sub("^.*_", "", expo_of_interest)
  expo_dat <- paste0(vac_interest, "_1")
  expo_level <- ifelse(dose_interest == "1", "2", "4")
  expo_mod <- paste0(expo_dat, expo_level)
  
  
  list(
    data = data_out,
    risk = risk,
    risk_null = risk_null,
    expo_dat = expo_dat,
    expo_level = expo_level,
    expo_mod = expo_mod,
    expo_interest = expo_of_interest,
    calendar_grp = calendar_adj
  )
}

## 1. Function to fit standard SCCS and fixed-effect Poisson model -------------
## for one exposure - outcome pair and one risk window

fit_sccs_3meth <- function(
    formula = event ~ bnt_1+ chad_1 + pre_bnt_1 + age, #The dependent variable should always be "event",
    indiv = id_num, 
    astart = obs_sta,
    aend = obs_end,
    aevent = diagnosis,
    adrug = list(cbind(bnt_1, bnt_2), 
                 cbind(chad_1, chad_2), 
                 cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)), 
    aedrug , 
    expogrp = list(c(0,1), # For BNT1, separate day 0 from risk window 
                   c(0,1), # For Chad1
                   c(0)),  # For pre-exposure window
    washout = list(), 
    sameexpopar = c(F, F, F), 
    agegrp, 
    seasongrp=NULL, dob=NULL, 
    dataformat="multi", 
    data,
    expo_of_interest,
    expo_mod,
    expo_dat,
    expo_level,
    outcome_name,
    model_type = "full_mod" # specify "null_mod" to fit a null model for Xu_2013
){
  # 1. Data reshaping (same as in function SCCS::standardsccs) -----------------
  yon <- deparse(substitute(adrug)) 
  yon1 <- as.formula(paste("z", "~", yon)) 
  adrugcolnames <- all.vars(yon1, functions = FALSE, unique = TRUE)[-1] # Get the names of all adrug columns
  
  adrug  <- eval(substitute(adrug), data, parent.frame())
  
  ## Changing adrug to a list if given as cbind(adrug1, adrug2,...) or adrug not as a list
  if ((dataformat=="multi" & !is.null(ncol(adrug)))) {
    adrug <- data.frame(adrug)
    adrug <- list(adrug) 
  } else if (dataformat=="stack" & !is.null(ncol(adrug))){
    adrug <- data.frame(adrug)
    adrug1 <- list()
    for (i in 1:ncol(adrug)){
      adrug1[[i]] <- adrug[,i]
    }
    adrug <- adrug1
  } else if (length(adrugcolnames)==1 & length(adrug)!=1) {
    adrug <- list(adrug)
    
  } else {
    adrug <- adrug
  }
  
  
  for (i in 1:length(adrug)){
    adrug[[i]] <- data.frame(adrug[[i]])
  }
  
  ncoladrug <- NULL
  for (i in 1:length(adrug)){
    ncoladrug[i] <- ncol(adrug[[i]]) 
  }
  
  
  for (i in 1:length(adrug)) {
    colnames(adrug[[i]]) <- adrugcolnames[c(1, cumsum(ncoladrug)+1)[-(length(ncoladrug)+1)][i]:cumsum(ncoladrug)[i]]
  }
  
  colname  <- adrugcolnames
  
  indiv  <- eval(substitute(indiv), data, parent.frame())
  astart <- eval(substitute(astart), data, parent.frame())
  aend   <- eval(substitute(aend), data, parent.frame())
  aevent <- eval(substitute(aevent), data, parent.frame())
  aedrug <- eval(substitute(aedrug), data, parent.frame())
  dob <- eval(substitute(dob), data, parent.frame())
  
  ## Changing aedrug to a list if given as cbind(aedrug1, aedrug2,...) or aedrug not as a list
  
  if ((dataformat=="multi" & !is.null(ncol(aedrug)))) {
    aedrug <- data.frame(aedrug)
    aedrug <- list(aedrug) 
  } else if (dataformat=="stack" & !is.null(ncol(aedrug))){
    aedrug <- data.frame(aedrug)
    aedrug1 <- list()
    for (i in 1:ncol(aedrug)){
      aedrug1[[i]] <- aedrug[,i]
    }
    aedrug <- aedrug1
  } else if (length(adrugcolnames)==1 & length(aedrug)!=1) {
    aedrug <- list(aedrug)
    
  } else {
    aedrug <- aedrug
  }
  
  
  ## Getting the fixed covariates from the formula  
  qq <- all.vars(as.formula(formula))[-c(which(all.vars(as.formula(formula))=="age"), which(all.vars(as.formula(formula))=="season"), which(all.vars(as.formula(formula))=="event"))]
  
  if (length(qq)==0) {
    cov <- cbind()
  }   else {
    cova <- qq[is.na(match(qq, colname))]
    cov <- data.frame(data[, cova])
    colnames(cov) <- cova
  }
  
  ## Reshape to SCCS data format
  
  chopdat <- SCCS::formatdata(indiv=indiv, astart=astart, aend=aend, aevent=aevent, 
                              adrug=adrug, aedrug=aedrug, expogrp = expogrp, washout = washout , 
                              sameexpopar = sameexpopar, agegrp = agegrp, seasongrp=seasongrp, 
                              dob=dob, cov=cov, dataformat=dataformat, data=NULL)
  
  ## Calculate the average time at risk
  exposure_data <- chopdat[chopdat[[expo_dat]] == expo_level, ]
  T_L <- sum(exposure_data$interval) / length(unique(exposure_data$indivL))
  T_L1 <- 1/ T_L
  
  
  # 2. Only fit fePois model if model_type = null_mod --------------------------
  
  # This is for the null model (to calculate the LR test statistic) for Xu_2013 
  
  if (model_type == "null_mod"){
    
    ## Make model formula
    fmla_null <- paste(formula, "+", "offset(log(interval))", "|", "indivL")
    fmla_null1 <- as.formula(paste("event~", fmla_null[3]))
    
    ## Fit the null fixed-effect Poisson regression and extract log-likelihood
    
    mod_null <- fixest::fepois(fml = fmla_null1, data = chopdat)
    sum_mod_null <- summary(mod_null)
    lr_null <- logLik(mod_null)
    
    
    model_stat_null <- data.frame(
      outcome = outcome_name,
      exposure = expo_of_interest,
      lr_null, 
      method = "fePois")
    
    return(model_stat_null)
    
  } else if (model_type == "full_mod") { 
    
  # 3. Fit standard SCCS model (for Xu_2011 and Campos_2017) -------------------
  
  ## Get the model formula and fit the model
  fmla <- paste(formula, "+", "strata(indivL)", "+", "offset(log(interval))")
  fmla1 <- as.formula(paste("event~", fmla[3]))
  mod <- clogit(formula = fmla1, data = chopdat)
  sum_mod <- summary(mod)
  
  ## Export model statistics
  main_expo <- expo_mod
  est_L <- stats::coef(sum_mod)[main_expo, 'coef'] 
  IRR_L <- stats::coef(sum_mod)[main_expo, 'exp(coef)']
  se_L <- stats::coef(sum_mod)[main_expo, 'se(coef)']
  p_val <- stats::coef(sum_mod)[main_expo, 'Pr(>|z|)']
  IRR_L_low_CI <- sum_mod$conf.int[main_expo, 'lower .95']
  IRR_L_up_CI <- sum_mod$conf.int[main_expo, 'upper .95']
  converge <- !(se_L > 10 | is.na(est_L) | is.nan(est_L) | is.infinite(se_L)) # Check convergence
  
  
  model_stat <- data.frame(
    outcome = outcome_name,
    exposure = expo_of_interest,
    T_L, T_L1, est_L, IRR_L,
    se_L, p_val,
    IRR_L_low_CI, IRR_L_up_CI,
    lr_full = NA, 
    method = "stdSCCS", converge)
  
  # 4. Fit fixed-effect Poisson models (for Xu_2013) ---------------------------
  
  ## Make model formula for the full model
  fmla_full <- paste(formula, "+", "offset(log(interval))", "|", "indivL")
  fmla_full1 <- as.formula(paste("event~", fmla_full[3]))
  
  
  ## Fit the full fixed-effect Poisson regression and extract log-likelihood
  
  mod_full <- fixest::fepois(fml = fmla_full1, data = chopdat)
  sum_mod_full <- summary(mod_full)
  lr_full <- logLik(mod_full)
  
  ## Export model statistics
  main_expo <- expo_mod
  est_L_fe <- sum_mod_full$coeftable[main_expo, 'Estimate'] 
  IRR_L_fe <- exp(est_L_fe)
  se_L_fe <- sum_mod_full$coeftable[main_expo, 'Std. Error']
  p_val_fe <- sum_mod_full$coeftable[main_expo, 'Pr(>|z|)']
  IRR_L_low_CI_fe <- exp(est_L_fe -1.96*se_L_fe)
  IRR_L_up_CI_fe <- exp(est_L_fe + 1.96*se_L_fe)
  converge_fe <- !(se_L_fe > 10 | is.na(est_L_fe) | is.nan(est_L) | is.infinite(se_L_fe))
  
  model_stat_fe <- data.frame(
    outcome = outcome_name,
    exposure = expo_of_interest,
    T_L = NA, T_L1 = NA,
    est_L = est_L_fe, IRR_L = IRR_L_fe,
    se_L = se_L_fe, p_val = p_val_fe,  
    IRR_L_low_CI = IRR_L_low_CI_fe, IRR_L_up_CI = IRR_L_up_CI_fe,
    lr_full, 
    method = "fePois", converge = converge_fe)
  
  
  all_results <- rbind(model_stat, model_stat_fe)
  
  return(all_results)
  }
  else {stop("Model type must be full_mod or null_mod")}
}


## 2. Function to loop `fit_sccs_3meth` through a series of candidate risk windows for one exposure -----

loop_risk_win <- function(
    data,
    risk_win,
    pre_ex = 30,
    max_risk_win = 70,
    expo,
    calendar_adjustment, 
    output_dir
){
  
  outcome_name <- tolower(unique(data$event_id))
  
  
  # Create the directory to store the results
  invisible(lapply(outcome_name, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  file_path <- file.path(output_dir, outcome_name, paste0(expo, ".csv"))
  
  # Loop through all candidate risk windows
  loop_risk_win_res <- foreach(k = seq_along(risk_win),
                               .combine = rbind) %do% 
    { 
      i <- risk_win[k]
      
      message(sprintf(
        "Scanning risk window: %d days (%s)",
        i, Sys.time()
      ))
      
      stage <- "sccs_prep"
      outcome_name <- NA_character_
      
      tryCatch({
        
        ## Stage 1: data preparation ----
        
        dat_sccs_prep <- sccs_prep(
          data = data,
          pre_ex = pre_ex,
          max_risk_win = max_risk_win,
          risk_scan = i,
          expo_of_interest = expo,
          calendar_adj = calendar_adjustment
        )
        
        outcome_name <- tolower(unique(dat_sccs_prep$data$event_id))
        
        ## Stage 2: model fitting ----
        
        ## 2.1.  Fit the standard SCCS model for Xu_2011 and Campos, and the full fePois model for Xu_2013
        stage <- "fit_sccs_3meth_fullmod"
        
        res_full <- fit_sccs_3meth(
          formula = event ~ bnt_1 + chad_1 + pre_bnt_1 + age,
          indiv = id_num,
          astart = obs_sta,
          aend = obs_end,
          aevent = diagnosis,
          adrug = list(
            cbind(bnt_1, bnt_2),
            cbind(chad_1, chad_2),
            cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)),
          aedrug = list(
            cbind(bnt_1 + dat_sccs_prep$risk[["bnt_1"]],
                  bnt_2 + dat_sccs_prep$risk[["bnt_2"]]),
            cbind(chad_1 + dat_sccs_prep$risk[["chad_1"]],
                  chad_2 + dat_sccs_prep$risk[["chad_2"]]),
            cbind(bnt_1 - 1, bnt_2 - 1, chad_1 - 1, chad_2 - 1)),
          expogrp = list(c(0,1), c(0,1), c(0)),
          washout = list(),
          sameexpopar = c(FALSE, FALSE, FALSE),
          agegrp = dat_sccs_prep$calendar_grp,
          data = dat_sccs_prep$data,
          expo_of_interest = dat_sccs_prep$expo_interest,
          expo_mod = dat_sccs_prep$expo_mod,
          expo_dat = dat_sccs_prep$expo_dat,
          expo_level = dat_sccs_prep$expo_level,
          outcome_name = outcome_name
        )
        
        ## Stage 2.2. Fit null model for Xu_2013
        
        stage <- "fit_sccs_3meth_nullmod"
        
        res_null <- fit_sccs_3meth(
          formula = event ~ bnt_1 + chad_1 + pre_bnt_1 + age,
          indiv = id_num,
          astart = obs_sta,
          aend = obs_end,
          aevent = diagnosis,
          adrug = list(
            cbind(bnt_1, bnt_2),
            cbind(chad_1, chad_2),
            cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)),
          aedrug = list(
            cbind(bnt_1 + dat_sccs_prep$risk_null[["bnt_1"]],
                  bnt_2 + dat_sccs_prep$risk_null[["bnt_2"]]),
            cbind(chad_1 + dat_sccs_prep$risk_null[["chad_1"]],
                  chad_2 + dat_sccs_prep$risk_null[["chad_2"]]),
            cbind(bnt_1 - 1, bnt_2 - 1, chad_1 - 1, chad_2 - 1)),
          expogrp = list(c(0,1), c(0,1), c(0)),
          washout = list(),
          sameexpopar = c(FALSE, FALSE, FALSE),
          agegrp = dat_sccs_prep$calendar_grp,
          data = dat_sccs_prep$data,
          expo_of_interest = dat_sccs_prep$expo_interest,
          expo_mod = dat_sccs_prep$expo_mod,
          expo_dat = dat_sccs_prep$expo_dat,
          expo_level = dat_sccs_prep$expo_level,
          outcome_name = outcome_name,
          model_type = "null_mod"
        )
        
        ## Stage 3. Merge results and compute LR test
        
        stage <- "merge_full_null"
        
        res_full2 <- res_full %>% 
          left_join(res_null, by = c("outcome", "exposure", "method")) %>%
          mutate(lr_test = lr_full - lr_null) %>%
          relocate(c(lr_null, lr_test), .after = lr_full)

        res_full2$candidate_risk_win <- i
        res_full2$max_risk_win <- max_risk_win
        res_full2$calendar_grp <- paste(calendar_adjustment, collapse = " ")
        
        # Render the results to a .csv file
        append_to_csv(res_full2, file_path)
        
        res_full2
        
      }, error = function(e) {
        
        ## Log error if occur ----
        log_error(
          sprintf(
            "Stage: %s | Exposure: %s | Risk window: %d | Outcome: %s | Error: %s",
            stage, expo, i, outcome_name, conditionMessage(e)
          )
        )
        
        ## ---- Placeholder result ----
        data.frame(
          outcome = outcome_name,
          exposure = expo,
          T_L = NA, T_L1 = NA,
          est_L = NA, IRR_L = NA,
          se_L = NA, p_val = NA,
          IRR_L_low_CI = NA, IRR_L_up_CI = NA,
          lr_full = NA, lr_null = NA, lr_test = NA,
          method = NA,
          converge = FALSE, 
          candidate_risk_win = NA,
          max_risk_win = max_risk_win, 
          calendar_grp = paste(calendar_adjustment, collapse = " ")
        )
      })
    }
}

  
## 3. Function to apply `loop_risk_win` through all exposures in a dataset -----

loop_4_exp <- function(
    data,
    risk_win = seq(1, 70, by = 1),
    pre_ex = 30,
    max_risk_win = 70,
    expo_to_scan = c("bnt_1", "bnt_2", "chad_1", "chad_2"),
    calendar_interval = 30,
    output_dir = here("Report", "Raw_results")
) {
  
  # Make sure that the intervals for the calendar time adjustment are within the observation period
  min_obs_sta <- min(data$obs_sta)
  max_obs_end <- max(data$obs_end)
  calendar <- seq(31, max_obs_end - (calendar_interval -1), by = calendar_interval)
  calendar_adj2 <- calendar[calendar > min_obs_sta]
  
  # Loop through all exposures
  
  foreach(expo = expo_to_scan, .combine = rbind) %do% {
    
    message(paste("Scanning for exposure:", expo, "-----------------------------"))
    
    res_1_expo <- loop_risk_win(data = data,
                                risk_win = risk_win,
                                pre_ex = pre_ex,
                                max_risk_win = max_risk_win,
                                expo = expo,
                                calendar_adjustment = calendar_adj2, 
                                output_dir = output_dir)
    res_1_expo
  }
  
}
  

################################################################################
# Part 2: Specific functions for Xu_2013 ---------------------------------------
################################################################################

## 2.1. Function to extract the optimal risk window with the maximum LRT -------
xu2013_maxlr <- function(
    result_table = results_raw_3meth  
){
  maxlr <- result_table %>% 
    filter(method == "fePois", converge == TRUE) %>%
    group_by(outcome, exposure) %>% 
    filter(lr_test == max(lr_test))
  
  return(maxlr)
}

## 2.2. Function to simulate data under the null hypothesis --------------------

sim_null_data_xu2013 <- function(
    seed,
    indiv = id_num, 
    astart = obs_sta, 
    aend = obs_end, 
    agegrp, 
    data,
    expo_of_interest
){
  
  set.seed(seed)
  
  # Step 1: Construct age groups & estimate dose-1 age distribution
  
  ## 1.1. Get information from data
  indiv_string <- deparse(substitute(indiv))
  astart_vec <- eval(substitute(astart), data, parent.frame())
  aend_vec   <- eval(substitute(aend), data, parent.frame())
  
  ## 1.2. Construct age-group intervals
  lower <- c(min(astart_vec), agegrp)
  upper <- c(agegrp - 1, max(aend_vec))
  
  agegrp_info <- data.table(
    lower = lower,
    upper = upper,
    age   = as.character(seq_along(lower)))
  
  ## 1.3. Extract vaccine exposure history
  
  vac_interest <- sub("_.*$", "", expo_of_interest)
  dose1 <- paste0(vac_interest,"_1")
  dose2 <- paste0(vac_interest,"_2")
  
  vacc_hist <- data[ , c(indiv_string, dose1, dose2)]
  
  names(vacc_hist) <- c("indiv", "vax_1", "vax_2")
  
  vacc_hist <- vacc_hist[!duplicated(vacc_hist$indiv), ]
  
  # Vaccination indicators
  vacc_hist$has_vax_1 <- !is.na(vacc_hist$vax_1)
  vacc_hist$has_vax_2 <- !is.na(vacc_hist$vax_2)
  
  # Inter-dose spacing (fixed)
  vacc_hist$delta12 <- with(
    vacc_hist,
    ifelse(has_vax_2, vax_2 - vax_1, NA))
  
  # Simulate dose 1 if the exposure of interest is dose 1: 
  
  if (expo_of_interest %in% c("bnt_1", "chad_1")){
    
    ## 1.4. Identify the age group of the dose 1 vaccination date
    # For each vaccination record, find the age-group interval 
    # whose [lower, upper] contains the vaccination age,
    #and assign that age-group number to the record. 
    
    vax1_obs <- vacc_hist[vacc_hist$has_vax_1, ]
    
    setDT(vax1_obs) # convert to data.table object
    
    vax1_obs[
      agegrp_info,
      age_group := i.age,
      on = .(vax_1 >= lower, vax_1 <= upper)]
    
    ## 1.5. Empirical age-group probabilities for dose 1
    prob_vax1_age <- prop.table(table(vax1_obs$age_group))
    
    
    # Step 2: randomly simulate dose 1 vaccination date for each individual using the above multinomial vaccination probabilities & shift dose 2
    
    
    vacc_hist$vax1_sim <- NA_integer_
    vacc_hist$vax2_sim <- NA_integer_
    
    for (i in which(vacc_hist$has_vax_1)) {

        ## 2.1. Randomly assign 'age' group at vaccination for dose 1
        age_grp <- sample(
          x    = names(prob_vax1_age),
          size = 1,
          prob = as.numeric(prob_vax1_age))
        
        bounds <- agegrp_info[age == age_grp]
        
        ## 2.2.Sample dose 1
        t1 <- round(runif(1, bounds$lower, bounds$upper))
        vacc_hist$vax1_sim[i] <- t1
        
        ## 2.3. Shift dose 2 date accordingly if present
        if (vacc_hist$has_vax_2[i]) {
          vacc_hist$vax2_sim[i] <- t1 + vacc_hist$delta12[i]
        }
        
        }
    
    # Step 3: format data to SCCS-compatible format
    # 3.1. Merge the simulated vaccination date to the original data 
    
    data <- merge(
      data,
      vacc_hist[, c("indiv", "vax1_sim", "vax2_sim")],
      by.x = indiv_string,
      by.y = "indiv",
      all.x = TRUE
    )
    
    # Only individuals with simulated values are overwritten
    if (expo_of_interest == "bnt_1") {
      data <- data %>%
        mutate(
          bnt_1 = ifelse(!is.na(vax1_sim), vax1_sim, bnt_1),
          bnt_2 = ifelse(!is.na(vax2_sim), vax2_sim, bnt_2))
    }
    
    if (expo_of_interest == "chad_1") {
      
      data <- data %>%
        mutate(
          chad_1 = ifelse(!is.na(vax1_sim), vax1_sim, chad_1),
          chad_2 = ifelse(!is.na(vax2_sim), vax2_sim, chad_2))
    }
    
  } else if (expo_of_interest %in% c("bnt_2", "chad_2")) {
    
    ## 1.4. Define weekly dosing interval bins and calculate the probability of each bin
    
    delta12_obs <- with(
      vacc_hist[vacc_hist$has_vax_2, ],
      vax_2 - vax_1
    )
    
    if (length(delta12_obs) == 0) {
      stop("No individuals with dose 2; cannot generate null distribution for dose 2")
    }
    
    delay_breaks <- seq(
      from = floor(min(delta12_obs)),
      to   = ceiling(max(delta12_obs)) + 7,
      by   = 7
    )
    
    # Assign the weekly bin to the observed dosing interval (e.g delay of 35 days belongs to bin [31, 38))
    delay_bins <- cut(
      delta12_obs,
      breaks = delay_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
    
    # multinomial probabilities for delay bins
    prob_delay_bin <- prop.table(table(delay_bins))
    
    # Step 2: randomly simulate dose 2 date using the above multinomial vaccination probabilities
    
    vacc_hist$vax2_sim <- NA_integer_
    
    for (i in which(vacc_hist$has_vax_2)) {
      
        
        ## 2.1. Randomly assign dosing interval 
        
        # sample delay bin
        bin_star <- sample(
          x    = names(prob_delay_bin),
          size = 1,
          prob = as.numeric(prob_delay_bin)
        )
        
        # extract numeric bounds of the bin
        bounds <- gsub("\\[|\\)|\\]", "", bin_star)
        bounds <- as.numeric(strsplit(bounds, ",")[[1]])
        
        # sample delay within the bin
        delta_star <- round(runif(1, min = bounds[1], max = bounds[2]))
        
        # Calculate new dose 2 date: = dose 1 + new interval
        vacc_hist$vax2_sim[i] <- vacc_hist$vax_1[i] + delta_star
        }
    
    # Step 3: format data to SCCS-compatible format
    # 3.1. Merge the simulated vaccination date to the original data 
    
    data <- merge(
      data,
      vacc_hist[, c("indiv", "vax2_sim")],
      by.x = indiv_string,
      by.y = "indiv",
      all.x = TRUE
    )
    
    if (expo_of_interest == "bnt_2") {
      data <- data %>%
        mutate(
          bnt_2 = ifelse(!is.na(vax2_sim), vax2_sim, bnt_2))
    }
    
    if (expo_of_interest == "chad_2") {
      data <- data %>%
        mutate(
          chad_2 = ifelse(!is.na(vax2_sim), vax2_sim, chad_2))
    }
    
  }
  return(data)
}

## 2.3. Function to get independent seeds for [10000] simulated datasets -------
get_seeds <- function(n_sim = 10000){
  set.seed(20251228, kind = "L'Ecuyer-CMRG", sample.kind = "Rejection")
  n_seed <- n_sim # Independent seed for each run in each scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}

## 2.4 Function to simulate the dataset and fit fePois model a large number of time ----

sim_null_dist_xu2013 <- function(seed_list,
                                 n_sim,
                                 data,
                                 expo_of_interest,
                                 calendar_adjustment,
                                 pre_ex = 30,
                                 max_risk_win = 70,
                                 result_table,
                                 output_dir){
  
  outcome_name <- tolower(unique(data$event_id))
  
  # Create the directory to store the results
  invisible(lapply(outcome_name, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  file_path <- file.path(output_dir, outcome_name, paste0(expo_of_interest, ".csv"))
  
  stage <- "xu2013_maxlr"
  optimal_win <- NA_real_
  
  ## Stage 1: get the optimal risk window for the exposure ----
  
  optimal_win <- xu2013_maxlr(result_table = result_table) %>%
    filter(
      exposure == expo_of_interest,
      outcome == outcome_name) %>%
    pull(candidate_risk_win)
  
  ## If there is no optimal_win, stop the procedure
  
  if (length(optimal_win) == 0 || is.na(optimal_win)) {
    message(
      sprintf(
        "Stage: xu2013_maxlr | Exposure: %s | Outcome: %s | No converged risk window",
        expo_of_interest, outcome_name))
    
    return(data.frame(
      outcome = outcome_name,
      exposure = expo_of_interest,
      optimal_win = NA,
      lr_full = NA_real_,
      lr_null = NA_real_,
      lr_test = NA_real_,
      method = "fePois", rep = NA, seed = NA))
  }
  
  
  null_dist <- foreach(i = 1:n_sim,
                       .combine = rbind) %do% 
    {
      if (i == 1 || i %% 20 == 0) {
        message(sprintf(
          "[%s] Replicate %d to %d out of %d",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          i, i + 19, n_sim))
      }
      
      tryCatch({
        
        ## Stage 2: simulate data under the null hypothesis for the respective exposure and outcome ----
        
        stage <- "sim_null_data_xu2013"
        
        sim_dat <- sim_null_data_xu2013(
          seed = seed_list[i],
          indiv = id_num,
          astart = obs_sta,
          aend = obs_end,
          agegrp = calendar_adjustment,
          data = data,
          expo_of_interest = expo_of_interest
        )
        
        ## Stage 3: prepare the data ----
        
        stage <- "sccs_prep"
        
        dat_sccs_prep <- sccs_prep(
          data = sim_dat,
          pre_ex = pre_ex,
          max_risk_win = max_risk_win,
          risk_scan = optimal_win,
          expo_of_interest = expo_of_interest,
          calendar_adj = calendar_adjustment
        )
        
        ## Stage 4: model fitting ----
        stage <- "fit_fePois_full"
        
        res_full <- fit_sccs_3meth(
          formula = event ~ bnt_1 + chad_1 + pre_bnt_1 + age,
          indiv = id_num,
          astart = obs_sta,
          aend = obs_end,
          aevent = diagnosis,
          adrug = list(
            cbind(bnt_1, bnt_2),
            cbind(chad_1, chad_2),
            cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)),
          aedrug = list(
            cbind(bnt_1 + dat_sccs_prep$risk[["bnt_1"]],
                  bnt_2 + dat_sccs_prep$risk[["bnt_2"]]),
            cbind(chad_1 + dat_sccs_prep$risk[["chad_1"]],
                  chad_2 + dat_sccs_prep$risk[["chad_2"]]),
            cbind(bnt_1 - 1, bnt_2 - 1, chad_1 - 1, chad_2 - 1)),
          expogrp = list(c(0,1), c(0,1), c(0)),
          washout = list(),
          sameexpopar = c(FALSE, FALSE, FALSE),
          agegrp = dat_sccs_prep$calendar_grp,
          data = dat_sccs_prep$data,
          expo_of_interest = dat_sccs_prep$expo_interest,
          expo_mod = dat_sccs_prep$expo_mod,
          expo_dat = dat_sccs_prep$expo_dat,
          expo_level = dat_sccs_prep$expo_level,
          outcome_name = outcome_name,
          model_type = "null_mod"
        )
        
        # Because this is the result for the full model, rename the column accordingly
        res_full <- res_full %>%
          rename(lr_full = lr_null)
        
        ## Stage 4.2. Fit null model
        
        stage <- "fit_fePois_nullmod"
        
        res_null <- fit_sccs_3meth(
          formula = event ~ bnt_1 + chad_1 + pre_bnt_1 + age,
          indiv = id_num,
          astart = obs_sta,
          aend = obs_end,
          aevent = diagnosis,
          adrug = list(
            cbind(bnt_1, bnt_2),
            cbind(chad_1, chad_2),
            cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)),
          aedrug = list(
            cbind(bnt_1 + dat_sccs_prep$risk_null[["bnt_1"]],
                  bnt_2 + dat_sccs_prep$risk_null[["bnt_2"]]),
            cbind(chad_1 + dat_sccs_prep$risk_null[["chad_1"]],
                  chad_2 + dat_sccs_prep$risk_null[["chad_2"]]),
            cbind(bnt_1 - 1, bnt_2 - 1, chad_1 - 1, chad_2 - 1)),
          expogrp = list(c(0,1), c(0,1), c(0)),
          washout = list(),
          sameexpopar = c(FALSE, FALSE, FALSE),
          agegrp = dat_sccs_prep$calendar_grp,
          data = dat_sccs_prep$data,
          expo_of_interest = dat_sccs_prep$expo_interest,
          expo_mod = dat_sccs_prep$expo_mod,
          expo_dat = dat_sccs_prep$expo_dat,
          expo_level = dat_sccs_prep$expo_level,
          outcome_name = outcome_name,
          model_type = "null_mod"
        )
        
        ## Stage 5. Merge results and compute LR test statistic
        
        stage <- "merge_full_null"
        
        res_full2 <- res_full %>% 
          left_join(res_null, by = c("outcome", "exposure", "method")) %>%
          mutate(lr_test = lr_full - lr_null, optimal_win = optimal_win) %>%
          relocate(c(lr_null, lr_test), .after = lr_full) %>%
          relocate(optimal_win, .after = exposure)
        
        res_full2$rep <- i
        res_full2$seed <- seed_list[i]
        
        # Export each simulated LR test stat to the .csv file
        append_to_csv(res_full2, file_path)
        res_full2
        
      }, error = function(e) {
        
        ## Log error if occur ----
        log_error(
          sprintf(
            "Stage: %s | Exposure: %s | Risk window: %d | Outcome: %s | Error: %s",
            stage, expo_of_interest, optimal_win, outcome_name, conditionMessage(e)
          )
        )
        
        ## ---- Placeholder result ----
        data.frame(
          outcome = outcome_name,
          exposure = expo_of_interest,
          optimal_win = NA,
          lr_full = NA, lr_null = NA, lr_test = NA,
          method = NA, rep = i, seed = seed_list[i]
        )
      })
      
    }
  return(null_dist)
  
}

## 2.5. Function to loop the simulation through 4 exposures ----

loop_4_exp_sim <- function(
    seed_list,
    n_sim,
    data,
    expo_to_scan = c("bnt_1", "bnt_2", "chad_1", "chad_2"),
    calendar_interval = 30,
    pre_ex = 30,
    max_risk_win = 70, 
    result_table, 
    output_dir = here("Report", "Null_dist")
) {

  
  # Make sure that the intervals for the calendar time adjustment are within the observation period
  min_obs_sta <- min(data$obs_sta)
  max_obs_end <- max(data$obs_end)
  calendar <- seq(31, max_obs_end - (calendar_interval -1), by = calendar_interval)
  calendar_adj2 <- calendar[calendar > min_obs_sta]
  
  foreach(expo = expo_to_scan, .combine = rbind) %do% {
    
    message(paste(sprintf("[%s] Simulate data for exposure",
                          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                  expo))
    
    sim_null_dist_xu2013(
      seed_list,
      n_sim,
      data,
      expo_of_interest = expo,
      calendar_adjustment = calendar_adj2,
      pre_ex = 30,
      max_risk_win = 70,
      output_dir = output_dir,
      result_table = result_table
    )
    
  }
  
}

## 2.6. Function to calculate the p-value for testing the null hypothesis 
## that there does not exist an interval with elevated risk

# For each exposure-outcome pair, this function extract the observed LR test stat, 
# load the respective simulated null distribution of the LR test, 
# calculate p-value and export the result to a .csv file

p_val_cal_xu2013 <- function(
    maxlr = xu2013_maxlr(result_table = results_raw_3meth),
    results_dir = here("Report", "Null_dist"),
    summary_dir = here("Report", "Summary"),
    p_val_file_name = "pval_Xu2013"
) {
  
  all_pval <- foreach(i = seq_len(nrow(maxlr)), .combine = rbind) %do% {
    
    outcome_i <- maxlr$outcome[i]
    expo_i    <- maxlr$exposure[i]
    lrt_obs_i <- maxlr$lr_test[i]
    
    sim_file <- paste0(expo_i, ".csv")
    file_path <- file.path(results_dir, outcome_i, sim_file)
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      return(NULL)
    }
    
    lrt_null <- tryCatch(
      read.csv(file_path),
      error = function(e) {
        warning(
          sprintf(
            "Failed to read null file | outcome=%s | exposure=%s | %s",
            outcome_i, expo_i, e$message
          ))
        return(NULL)
      })
    
    if (is.null(lrt_null) || nrow(lrt_null) == 0) {
      return(NULL)
    }
    
    pval_i <- tryCatch({
      B <- nrow(lrt_null)
      (1 + sum(lrt_null$lr_test_sim >= lrt_obs_i)) / (1 + B)
    }, error = function(e) {
      warning(
        sprintf(
          "Error calculating p-value | outcome=%s | exposure=%s | %s",
          outcome_i, expo_i, e$message
        ))
      return(NULL)
    })
    
    if (is.null(pval_i)) {
      return(NULL)
    }
    
    data.frame(
      outcome  = outcome_i,
      exposure = expo_i,
      p_value  = pval_i
    )
  }
  
  output_path <- file.path(
    summary_dir,
    paste0(p_val_file_name, ".csv")
  )
  
  readr::write_csv(all_pval, output_path)
  
  all_pval
}


################################################################################
# Part 3: Functions to select the optimal window per method --------------------
################################################################################

## 3.1. Function to select the optimal window for Xu_2011 (per exposure-outcome pair) ----
# This function automatically selects window with max IRR and produce the plot of IRR_L vs 1/T(L)
# It is necessary to check if the plot shows an approximate linear relationship 
# between IRRL and 1/T(L) when L>LM

Xu_2011_plot <- function(data,
                         x = "T_L1",
                         y = "IRR_L",
                         xlab = "1/ T(L)",
                         ylab = "R (L)", 
                         xbreak = 20,
                         plot_dir = file.path(here("Report"), "Plot")
                         ) {
  
  ## Early exit if data is empty
  if (is.null(data) || nrow(data) == 0) {
    return(
      tibble::tibble(
        outcome = character(),
        exposure = character(),
        est_L = numeric(),
        se_L = numeric(),
        IRR_L = numeric(),
        IRR_L_low_CI = numeric(),
        IRR_L_up_CI = numeric(),
        p_val = numeric(),
        lr_test = numeric(),
        optimal_window_days = numeric(),
        method = character(),
        error = character()
        
      )
    )
  }
  
  
  create_directory(plot_dir)
  
  # Only select candidate risk windows with converged SCCS model
  data2 <- data %>% 
    filter(method == "stdSCCS",
           converge == TRUE)
  
  if (nrow(data2) < 2) {
    return(
      tibble::tibble(
        outcome = character(),
        exposure = character(),
        est_L = numeric(),
        se_L = numeric(),
        IRR_L = numeric(),
        IRR_L_low_CI = numeric(),
        IRR_L_up_CI = numeric(),
        p_val = numeric(),
        lr_test = numeric(),
        optimal_window_days = numeric(),
        method = character(),
        error = "insufficient data"
      )
    )
  }
  
  # row with maximum IRR
  max_row <- data2 %>% 
    slice_max(.data[[y]], n = 1, with_ties = FALSE)
  
  # Make plot of IRR_L versus 1/T(L) (average time at risk)
  
  plot_name <- paste0("Xu2011_", max_row$outcome, "_", max_row$exposure, ".png")
  
  png(file.path(plot_dir, plot_name), 
      width = 20, height = 12, units = "cm", res = 300)
  
  plot_irr <- ggplot(data2, aes(x = .data[[x]], y = .data[[y]])) +
    geom_line() +
    geom_point() +
    geom_point(
      data = max_row,
      size = 2,
      colour = "darkred"
    ) +
    scale_x_continuous(
      name = xlab,
      breaks = scales::breaks_pretty(n = xbreak),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      name = ylab,
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(subtitle = paste0("Max IRR = ", round(max_row$IRR_L, 2),  " at L = ", max_row$candidate_risk_win, " days")) + 
    theme_minimal()
  
  print(plot_irr)
  dev.off()
  
  res_xu2011 <- max_row %>%
    select(outcome, exposure, est_L, se_L, IRR_L, IRR_L_low_CI, IRR_L_up_CI, p_val) 
  
  res_xu2011 <- res_xu2011 %>%
    mutate(lr_test = NA,
      optimal_window_days = max_row$candidate_risk_win,
      method = "Xu_2011",
      error = NA_character_)
  
  return(res_xu2011)
}

## 3.2. Function to select the optimal window for Campos (per exposure-outcome pair) -----

# This Function makes a scatterplot/line plot of the estimated IRR and the 
# (inverse of) a sequence of risk period lengths, and fit a sequence of linear-quadratic 
# spline models for different knot points through the scatter plot
 
Campos_plot <- function(data, 
                        risk_win = "candidate_risk_win",# variable 'risk length' in the data
                        T_L_inv = "T_L1", # variable '1/risk length' in the data
                        IRR_L = "IRR_L",   # variable 'estimated IRR' in the data
                        xbreak = 20, # number of breaks for the x axis plot
                        plot_dir = file.path(here("Report"), "Plot")
){
  ## Early exit if data is empty
  if (is.null(data) || nrow(data) == 0) {
    return(
      tibble::tibble(
        outcome = character(),
        exposure = character(),
        est_L = numeric(),
        se_L = numeric(),
        IRR_L = numeric(),
        IRR_L_low_CI = numeric(),
        IRR_L_up_CI = numeric(),
        p_val = numeric(),
        lr_test = numeric(),
        optimal_window_days = numeric(),
        method = character(),
        error = character()
      )
    )
  }
  
  create_directory(plot_dir)
  
  # Only select candidate risk windows with converged SCCS model
  data <- data %>% 
    filter(method == "stdSCCS",
           converge == TRUE)
  
  if (nrow(data) < 3) { # avoid unstable spline fit
    return(
      tibble::tibble(
        outcome = data$outcome[1],
        exposure = data$exposure[1],
        est_L = NA_real_,
        se_L = NA_real_,
        IRR_L = NA_real_,
        IRR_L_low_CI = NA_real_,
        IRR_L_up_CI = NA_real_,
        p_val = NA_real_,
        lr_test = NA_real_,
        optimal_window_days = NA_real_,
        method = "Campos",
        error = "insufficient data"
      )
    )
  }
  
  # 1. Define range of candidate knots
  # Convert to the inverse scale (1/t) for the x-axis variable
  candidate_risk_win <- data[[risk_win]]
  knot_candidates_inv <- data[[T_L_inv]] #only use knots as the candidate risk windows in the SCCS
  
  # Get data for fitting splines
  x <- data[[T_L_inv]]
  y <- data[[IRR_L]]
  
  # 2. Fit linear-quadratic splines across candidate knots and extract sum of square error
  ## Define function  
  calc_sse <- function(knot_inv){ #knot_inv = 1/risk length
    
    # Term for the change in slope after the knot:
    lp  <- pmax(0, x - knot_inv) # = 0 before the knot
    # Term for the curvature after the knot
    lp2 <- lp^2
    # Fit linear-quadratic spline
    lqs <- lm(y ~ x + lp + lp2) 
    # Calculate SSE
    sum(resid(lqs)^2)
  }
  ## Apply function
  sse_values <- sapply(knot_candidates_inv, calc_sse)
  
  # 3. Find the knot (in inverse scale) that minimizes SSE
  best_knot_inv <- knot_candidates_inv[which.min(sse_values)]
  
  
  # 4. Fit final spline using the best knot
  lp  <- pmax(0, x - best_knot_inv)
  lp2 <- lp^2
  
  final_model <- lm(y ~ x + lp + lp2)
  
  # 5. Generate prediction grid for plotting
  x_grid <- seq(min(x), max(x), length.out = 200)
  
  grid_df <- data.frame(
    x   = x_grid,
    lp  = pmax(0, x_grid - best_knot_inv),
    lp2 = pmax(0, x_grid - best_knot_inv)^2
  )
  
  grid_df$pred_y <- predict(final_model, newdata = grid_df)
  
  # 6. Visualise the results
  optimal_point <- data[data[[T_L_inv]] == best_knot_inv, ]
  
  plot_name <- paste0("Campos_", optimal_point$outcome, "_", optimal_point$exposure, ".png")
  
  png(file.path(plot_dir, plot_name), 
      width = 20, height = 12, units = "cm", res = 300)
  
  p <- ggplot(data, aes(x = .data[[T_L_inv]], y = .data[[IRR_L]])) +
    geom_point(alpha = 0.8, color = "black") +
    geom_line(data = grid_df, aes(x = x, y = pred_y), linewidth = 1) +
    geom_point(data = optimal_point, 
               aes(x = .data[[T_L_inv]], .data[[IRR_L]]),
               color = "darkred", 
               size = 3) + 
    labs(title = "Optimal Quadratic Spline Fit (Campos 2017)",
         subtitle = paste("Optimal knot at", optimal_point$candidate_risk_win, "days"), 
         y = "IRR (L)") +
    scale_x_continuous(
      name = "Inverse Risk Length (1/T(L))",
      breaks = scales::breaks_pretty(n = xbreak)
    ) +
    theme_minimal()
  
  print(p)
  dev.off()
  
  # Export the result of the optimal window
  res_campos <- optimal_point %>%
    select(outcome, exposure, est_L, se_L, IRR_L, IRR_L_low_CI, IRR_L_up_CI, p_val, lr_test) 
  
  res_campos <- res_campos %>%
    mutate(
           optimal_window_days = optimal_point$candidate_risk_win,
           method = "Campos",
           error = NA_character_)
  
  return(res_campos)
}

## No function to select the optimal window for Xu_2013 because we just need to select the one with max LR test stat
## (will perform it directly in the loop)
    
## 3.3 Function to apply the selection of optimal window by all methods and for all exposure-outcome pairs ----

# This function will return the optimal risk windows selected by three methods for each exposure-outcome pair, 
# and export to a .csv file

select_risk_win <- function(data, # the table with results from all SCCS models with all candidate windows
                            summary_dir = here("Report", "Summary"),
                            plot_dir = here("Report", "Plot"),
                            summary_file_name = "Summary_all_scens") {
  
  create_directory(summary_dir)
  
  ## Split data
  data_split <- data %>%
    dplyr::group_by(outcome, exposure) %>%
    dplyr::group_split()
  
  all_res <- purrr::map_dfr(data_split, function(data2) {
    
    outcome_i  <- data2$outcome[1]
    exposure_i <- data2$exposure[1]
    
    ## Xu 2011 ----
    res_xu2011 <- tryCatch(
      Xu_2011_plot(data = data2, plot_dir = plot_dir),
      error = function(e) {
        log_error(
          sprintf(
            "Xu_2011 failed | outcome=%s | exposure=%s | %s",
            outcome_i, exposure_i, conditionMessage(e))
        )
        NULL
      })
    
    ## Campos ----
    res_campos <- tryCatch(
      Campos_plot(data = data2, plot_dir = plot_dir),
      error = function(e) {
        log_error(
          sprintf(
            "Campos failed | outcome=%s | exposure=%s | %s",
            outcome_i, exposure_i, conditionMessage(e)
          ))
        NULL
      })
    
    ## Xu 2013 ----
    res_xu2013 <- tryCatch({
      
      base <- data2 %>%
        dplyr::filter(method == "fePois", converge == TRUE)
      
      if (nrow(base) == 0) {
        log_error(
          sprintf(
            "Xu_2013 skipped (no valid rows) | outcome=%s | exposure=%s",
            outcome_i, exposure_i
          ))
        return(NULL)
      }
      
      if (nrow(base) < 2) {
        return(
          tibble::tibble(
            outcome = data$outcome[1],
            exposure = data$exposure[1],
            est_L = NA_real_,
            se_L = NA_real_,
            IRR_L = NA_real_,
            IRR_L_low_CI = NA_real_,
            IRR_L_up_CI = NA_real_,
            p_val = NA_real_,
            lr_test = NA_real_,
            optimal_window_days = NA_real_,
            method = "Xu_2013",
            error = "insufficient data"
          )
        )
      }
      
      best_lr <- max(base$lr_test, na.rm = TRUE)
      
      base %>%
        dplyr::filter(lr_test == best_lr) %>%
        dplyr::select(
          outcome, exposure,
          est_L, se_L,
          IRR_L, IRR_L_low_CI, IRR_L_up_CI,
          p_val, lr_test, candidate_risk_win
        ) %>%
        dplyr::rename(optimal_window_days = candidate_risk_win) %>%
        dplyr::mutate(method = "Xu_2013",
                      error = NA_character_)
      
    }, error = function(e) {
      log_error(
        sprintf(
          "Xu_2013 failed | outcome=%s | exposure=%s | %s",
          outcome_i, exposure_i, conditionMessage(e)
        ))
      NULL
    })
    
    ## Bind results
    dplyr::bind_rows(
      res_xu2011,
      res_campos,
      res_xu2013
    )
  })
  
  ## Export summary
  output_path <- file.path(
    summary_dir,
    paste0(summary_file_name, ".csv")
  )
  
  readr::write_csv(all_res, output_path)
  
  all_res
}



  

  
  


