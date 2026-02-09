

################################################################################
# Part 1: Define functions to implement three methods --------------------------
################################################################################

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

# Function to compute pre-exposure period for all exposures in the dataset:

sccs_prep <- function(data,
                      risk_scan, # The risk window we are scanning for the expo_of_interest
                      expo_of_interest, # The exposure we are scanning the risk window for,
                      pre_ex = 30, # Length of pre-exposure period
                      max_risk_win = 70, # The default risk window length for all exposures
                      calendar_adj
                       ) {
  
  # Assign correct risk window to each exposure
  risk <- c(
    scan   = risk_scan,
    bnt_1  = max_risk_win,
    bnt_2  = max_risk_win,
    chad_1 = max_risk_win,
    chad_2 = max_risk_win
  )
  
  risk[expo_of_interest] <- risk["scan"]
  
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
  
  # Obtain variable names for fitting model
  vac_interest <- sub("_.*$", "", expo_of_interest)
  dose_interest <- sub("^.*_", "", expo_of_interest)
  expo_dat <- paste0(vac_interest, "_1")
  expo_level <- ifelse(dose_interest == "1", "2", "4")
  expo_mod <- paste0(expo_dat, expo_level)
  
  
  list(
    data = data_out,
    risk = risk,
    expo_dat = expo_dat,
    expo_level = expo_level,
    expo_mod = expo_mod,
    expo_interst = expo_of_interest,
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
    outcome_name
){
  # 1. Data reshaping -----------------------------------------------------
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
  
  # 2. Fit standard SCCS model (for Xu_2011 and Campos_2017) -------------------
  
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
  converge <- !(se_L > 10 | is.na(se_L) | is.nan(se_L) | is.infinite(se_L)) # Check convergence
  
  
  model_stat <- data.frame(
    outcome = outcome_name,
    exposure = expo_of_interest,
    T_L, T_L1, est_L, IRR_L,
    se_L, p_val,
    IRR_L_low_CI, IRR_L_up_CI,
    lr_null = NA, lr_full = NA, lr_test = NA, 
    method = "stdSCCS", converge)
  
  # 3. Fit fixed-effect Poisson models (for Xu_2013) ---------------------------
  
  ## Make model formula for the full model
  fmla_full <- paste(formula, "+", "offset(log(interval))", "|", "indivL")
  fmla_full1 <- as.formula(paste("event~", fmla_full[3]))
  
  ## Make model formula for the null model (model without the effect that we want to scan for the risk window)
  #fmla_null <- paste(formula_null, "+", "offset(log(interval))", "|", "indivL")
  #fmla_null1 <- as.formula(paste("event~", fmla_null[3]))
  
  ## Fit the full and null fixed-effect Poisson regression and extract log-likelihood
  
  mod_full <- fixest::fepois(fml = fmla_full1, data = chopdat)
  sum_mod_full <- summary(mod_full)
  lr_full <- logLik(mod_full)
  
  #mod_null <- fixest::fepois(fml = fmla_null1, data = chopdat)
  #lr_null <- logLik(mod_null)
  
  ## log likelihood ratio test statistic
  #lr_test <- lr_full - lr_null
  
  ## Export model statistics
  main_expo <- expo_mod
  est_L_fe <- sum_mod_full$coeftable[main_expo, 'Estimate'] 
  IRR_L_fe <- exp(est_L_fe)
  se_L_fe <- sum_mod_full$coeftable[main_expo, 'Std. Error']
  p_val_fe <- sum_mod_full$coeftable[main_expo, 'Pr(>|z|)']
  IRR_L_low_CI_fe <- exp(est_L_fe -1.96*se_L_fe)
  IRR_L_up_CI_fe <- exp(est_L_fe + 1.96*se_L_fe)
  converge_fe <- !(se_L_fe > 10 | is.na(se_L_fe) | is.nan(se_L_fe) | is.infinite(se_L_fe))
  
  model_stat_fe <- data.frame(
    outcome = outcome_name,
    exposure = expo_of_interest,
    T_L = NA, T_L1 = NA,
    est_L = est_L_fe, IRR_L = IRR_L_fe,
    se_L = se_L_fe, p_val = p_val_fe,  
    IRR_L_low_CI = IRR_L_low_CI_fe, IRR_L_up_CI = IRR_L_up_CI_fe,
    lr_null = NA, lr_full , lr_test = NA, 
    method = "fePois", converge = converge_fe)
  
  
  all_results <- rbind(model_stat, model_stat_fe)
  
  return(all_results)
}

## 2. Function to loop `fit_sccs_3meth` through 4 exposures --------------------

loop_4_exp <- function(
    data,
    risk_scan,
    pre_ex = 30,
    max_risk_win = 70,
    expo_to_scan = c("bnt_1", "bnt_2", "chad_1", "chad_2"),
    calendar_adjustment
) {
  
  foreach(expo = expo_to_scan, .combine = rbind) %do% {
    
    dat_sccs_prep <- sccs_prep(
      data = data,
      pre_ex = pre_ex,
      max_risk_win = max_risk_win,
      risk_scan = risk_scan,
      expo_of_interest = expo,
      calendar_adj = calendar_adjustment
    )
    
    fit_sccs_3meth(
      formula = event ~ bnt_1 + chad_1 + pre_bnt_1 + age,
      indiv = id_num,
      astart = obs_sta,
      aend = obs_end,
      aevent = diagnosis,
      adrug = list(
        cbind(bnt_1, bnt_2),
        cbind(chad_1, chad_2),
        cbind(pre_bnt_1, pre_bnt_2, pre_chad_1, pre_chad_2)
      ),
      aedrug = list(
        cbind(bnt_1 + dat_sccs_prep$risk[["bnt_1"]],
              bnt_2 + dat_sccs_prep$risk[["bnt_2"]]),
        cbind(chad_1 + dat_sccs_prep$risk[["chad_1"]],
              chad_2 + dat_sccs_prep$risk[["chad_2"]]),
        cbind(bnt_1 - 1, bnt_2 - 1, chad_1 - 1, chad_2 - 1)
      ),
      expogrp = list(c(0,1), c(0,1), c(0)),
      washout = list(),
      sameexpopar = c(FALSE, FALSE, FALSE),
      agegrp = dat_sccs_prep$calendar_grp,
      data = dat_sccs_prep$data,
      expo_of_interest = dat_sccs_prep$expo_interst,
      expo_mod = dat_sccs_prep$expo_mod,
      expo_dat = dat_sccs_prep$expo_dat,
      expo_level = dat_sccs_prep$expo_level,
      outcome_name = tolower(unique(dat_sccs_prep$data$event_id))
    )
  }
}

  
## 3. Function to apply `loop_4_exp` through a series of risk windows in a dataset ----

loop_risk_win <- function(risk_win = seq(1, 70, by = 1), 
                          data,
                          pre_ex = 30,
                          calendar_interval = 30,
                          max_risk_win = 70
                          ){
  
  
  # Make sure that the intervals for the calendar time adjustment are within the observation period
  min_obs_sta <- min(data$obs_sta)
  max_obs_end <- max(data$obs_end)
  calendar <- seq(31, max_obs_end - (calendar_interval -1), by = calendar_interval)
  calendar_adj2 <- calendar[calendar > min_obs_sta]
  
  # Loop through all candidate risk windows
  loop_risk_win_res <- foreach(k = seq_along(risk_win),
                               .combine = rbind) %do% 
    { 
      i <- risk_win[k]
      message(paste("Scanning risk window:", i, "days, at", Sys.time()))
      
      res_1_risk_win <- loop_4_exp(risk_scan = i, 
                                   data = data, 
                                   max_risk_win = max_risk_win,
                                   calendar_adjustment = calendar_adj2)
      res_1_risk_win$candidate_risk_win <- i
      res_1_risk_win$max_risk_win <- max_risk_win
      res_1_risk_win$calendar_grp <- paste(calendar_adj2, collapse = " ")
      res_1_risk_win
      }
}



