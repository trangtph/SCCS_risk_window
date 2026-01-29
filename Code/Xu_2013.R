#############################################################
# Analysis code for Xu (2013) --------------- ###############
#############################################################

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)
pacman::p_load(
  SCCS,
  dplyr,
  survival,
  foreach,
  ggplot2,
  future,
  doFuture,
  tictoc,
  progressr
)


# The method of Xu (2013) comprises two steps:

# Step 1: Fit a series of fixed-effect Poisson models and compute the likelihood ratio test ----------------

## 1.1. Function Xu_2013 performs the following steps: ----
# (1) Data preprocessing and reshaping: same as in function `Xu_2011`
#     - Processes input variables and reformats the data into SCCS-compatible
#       interval format using SCCS::formatdata.
#     - This follows the same internal procedure as in `SCCS::standardsccs`.
# (2) Fit the full fixed-effect Poisson model with all vaccine and age effect 
#     The fixed-effect Poisson model (function `fepois`) is equivalent to the 
#     conditional Poisson model (Xu et al., J Biom Biostat 2012)
# (3) Fit the null Poisson model without the effect of the vaccine that we are scanning the risk window for
#     The likelihood of this model is the likelihood function under the null hypothesis 
#     that the incidence rate of adverse events remains constant throughout the follow-up 
#     period while adjusting for other covariates
# (4) Compute the likelihood ratio test (LRT) statistic: 
#     LRT = loglikelihood(full model) - loglikelihood(null model)
#     The risk window corresponding to the maximum LRT is the optimal risk window

Xu_2013 <- function(
    formula, #The dependent variable should always be "event",
    formula_null = NULL, # formula without the exposure of interest
    indiv, astart, aend, aevent, adrug, aedrug, expogrp = list(), washout = list(), 
    sameexpopar = list(), agegrp = NULL, seasongrp=NULL, dob=NULL, dataformat="stack", data,
    expo_of_interest, risk_win_of_interest = 1
) 
{
  # 1. Data reshaping -----------------------------------------------------
  yon <- deparse(substitute(adrug)) 
  yon1 <- as.formula(paste("z", "~", yon)) 
  adrugcolnames <- all.vars(yon1, functions = FALSE, unique = TRUE)[-1] # Get the names of all adrug columns
  
  adrug  <- eval(substitute(adrug), data, parent.frame())
  
  # Changing adrug to a list if given as cbind(adrug1, adrug2,...) or adrug not as a list
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
  
  # Changing aedrug to a list if given as cbind(aedrug1, aedrug2,...) or aedrug not as a list
  
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
  
  
  # Getting the fixed covariates from the formula
  
  qq <- all.vars(as.formula(formula))[-c(which(all.vars(as.formula(formula))=="age"), which(all.vars(as.formula(formula))=="season"), which(all.vars(as.formula(formula))=="event"))]
  
  if (length(qq)==0) {
    cov <- cbind()
  }   else {
    cova <- qq[is.na(match(qq, colname))]
    cov <- data.frame(data[, cova])
    colnames(cov) <- cova
  }
  
  # Reshape to SCCS data format
  
  chopdat <- SCCS::formatdata(indiv=indiv, astart=astart, aend=aend, aevent=aevent, 
                              adrug=adrug, aedrug=aedrug, expogrp = expogrp, washout = washout , 
                              sameexpopar = sameexpopar, agegrp = agegrp, seasongrp=seasongrp, 
                              dob=dob, cov=cov, dataformat=dataformat, data=NULL)

  # 2. Fit fixed-effect Poisson model ------------------------------------------
  # Make model formula for the full model
  fmla_full <- paste(formula, "+", "offset(log(interval))", "|", "indivL")
  fmla_full1 <- as.formula(paste("event~", fmla_full[3]))
  
  # Make model formula for the null model (model without the effect that we want to scan for the risk window)
  fmla_null <- paste(formula_null, "+", "offset(log(interval))", "|", "indivL")
  fmla_null1 <- as.formula(paste("event~", fmla_null[3]))
  
  # Fit the full and null fixed-effect Poisson regression and extract log-likelihood
  
  mod_full <- fixest::fepois(fml = fmla_full1, data = chopdat)
  sum_mod_full <- summary(mod_full)
  lr_full <- logLik(mod_full)
  
  mod_null <- fixest::fepois(fml = fmla_null1, data = chopdat)
  lr_null <- logLik(mod_null)
  
  # log likelihood ratio test statistic
  lr_test <- lr_full - lr_null
  
  # Export model statistics
  main_expo <- paste0(expo_of_interest, risk_win_of_interest)
  est_L <- sum_mod_full$coeftable[main_expo, 'Estimate'] 
  IRR_L <- exp(est_L)
  se_L <- sum_mod_full$coeftable[main_expo, 'Std. Error']
  IRR_L_low_CI <- exp(est_L -1.96*se_L)
  IRR_L_up_CI <- exp(est_L + 1.96*se_L)
  
  model_stat <- data.frame(est_L, IRR_L,
                           se_L, IRR_L_low_CI, IRR_L_up_CI,
                           lr_null, lr_full, lr_test)
  
  return(model_stat) 

}

## 1.2. Function to extract the  optimal risk window with the maximum LRT -----------
xu2013_maxlr <- function(
    result_table  
){
  maxlr <- result_table[result_table$lr_test == max(result_table$lr_test),]
  return(maxlr)
}

# Step 2: Test the null hypothesis ----
# that there does not exist a risk window after vaccination with elevated risk of adverse events

# We use Monte Carlo simulation to obtain the p-value: 
# - Simulate the vaccination date under the null hypothesis that there does not 
#   exist a risk window after vaccination with elevated risk of adverse events:
#   the vaccination date only depends on the age/ calendar time, and is not 
#   associate with the event date
# - Use the observed event date in our original data and the simulated vaccination date 
#   as the new dataset, reshape this to SCCS-compatible format
# - Fit the full and the null model (similar to those in the function `Xu_2013` above) to the new dataset
# - The risk window in the full model is the optimal risk window identified from Step 1
# - Obtain the LRT
# - Repeat the above four steps a large number of time (10000 times), we obtain the
#   distribution of the LRT statistic under the null hypothesis
# - The p-value is calculated as the number of time the null test statistics equal
#   or larger than the observed test statistic corresponding to the optimal risk window
#   (obtained from the function `xu2013_maxlr`)

## 2.1. Function to simulate the LRT statistic under the null hypothesis -------------
lr_test_sim_null <- function(
    rep = 1,
    seed,
    data,
    maxlr = xu2013_maxlr(),
    formula, #The dependent variable should always be "event", one of the exposure should always be "vax_date_sim"
    formula_null = NULL, # formula without the exposure of interest
    indiv, astart, aend, aevent, 
    adrug,  # one of the exposure should always be "vax_date_sim" 
    aedrug, # risk window of vax_date_sim should = the "optimal" risk window identify by the Xu_2013 method
    expogrp = list(), agegrp = NULL,
    expo_of_interest 
){
  # Step 1: obtain the vaccination probabilities across the 'age' groups
  # 1.1. Construct the table of the day interval for each 'age' group
  astart <- eval(substitute(astart), data, parent.frame())
  aend   <- eval(substitute(aend), data, parent.frame())
  lower <- c(min(astart), agegrp)
  upper <- c(agegrp-1, max(aend))
  
  agegrp_info <- data.table(lower, upper, age = as.character(seq_along(lower)))
  
  # 1.2. Get the vaccination date per individual
  indiv_string <- deparse(substitute(indiv))
  data_vax <- data[!duplicated(data[, c(indiv_string, expo_of_interest)]), 
                   c(indiv_string, expo_of_interest)]
  names(data_vax) <- c("indiv", "vax")
  
  setDT(data_vax) # convert to data.table object
  
  # 1.3. Identify the age group of the vaccination date
  # For each vaccination record, find the age-group interval 
  # whose [lower, upper] contains the vaccination age,
  #and assign that age-group number to the record.
  data_vax[
    agegrp_info,
    age_group := i.age,
    on = .(vax >= lower, vax <= upper)
  ]
  
  # 1.4. Calculate vaccination probability for each 'age' group
  prob_vax_age <- prop.table(table(data_vax$age_group))
  
  # Step 2: randomly simulate vaccination date for each individual using the above multinomial vaccination probabilities
  # 2.1.Randomly assign 'age' group at vaccination
  set.seed(seed)
  
  data_vax$vax_age_sim <- sample(x = names(prob_vax_age),
                                 size = nrow(data_vax),
                                 replace = TRUE,
                                 prob = as.numeric(prob_vax_age))
  # 2.2. Attach 'age' intervals
  
  data_vax <- merge(data_vax, as.data.frame(agegrp_info),
                    by.x = "vax_age_sim", by.y = "age", all.x = TRUE)
  
  # 2.3. Simulate vaccination date within 'age' interval
  data_vax$vax_date_sim <- with(data_vax, 
                                round(runif(nrow(data_vax), min = lower, max = upper)))
  
  data_vax <- data_vax %>% select(indiv, vax_date_sim)
  
  
  # Step 3: format data to SCCS-compatible format
  # 3.1. Merge the simulated vaccination date to the original data 
  
  data <- merge(data, data_vax, by.x = indiv_string, by.y = "indiv", all.x = TRUE)
  
  # 3.2. Reformat the original data 
  
  # "vax_date_sim" being the main exposure
  # risk window for "vax_date_sim" = the "optimal" risk window identify by the Xu_2013 method
  
  # This part is the same as in the function "Xu_2013"
  
  yon <- deparse(substitute(adrug)) 
  yon1 <- as.formula(paste("z", "~", yon)) 
  adrugcolnames <- all.vars(yon1, functions = FALSE, unique = TRUE)[-1] # Get the names of all adrug columns

  adrug  <- eval(substitute(adrug), data, parent.frame())
  
  dataformat <- "stack"
  # Changing adrug to a list if given as cbind(adrug1, adrug2,...) or adrug not as a list
  if (dataformat=="stack" & !is.null(ncol(adrug))){
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
  aevent <- eval(substitute(aevent), data, parent.frame())
  aedrug <- eval(substitute(aedrug), data, parent.frame()) 
  
  if (dataformat=="stack" & !is.null(ncol(aedrug))){
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
  
  # Reshape to SCCS data format
  
  chopdat <- SCCS::formatdata(indiv=indiv, astart=astart, aend=aend, aevent=aevent, 
                              adrug=adrug, aedrug=aedrug, expogrp = list(), washout = list() , 
                              sameexpopar = list(), agegrp = agegrp, seasongrp=NULL, 
                              dob=NULL, cov=cbind(), dataformat="stack", data=NULL)
  
  # Step 4: Fit full and null model and extract the log-likelihood ratio test
  # Make model formula for the full model
  fmla_full <- paste(formula, "+", "offset(log(interval))", "|", "indivL")
  fmla_full1 <- as.formula(paste("event~", fmla_full[3]))
  
  # Make model formula for the null model (model without the effect that we want to scan for the risk window)
  fmla_null <- paste(formula_null, "+", "offset(log(interval))", "|", "indivL")
  fmla_null1 <- as.formula(paste("event~", fmla_null[3]))
  
  # 4.1. Fit full model
  mod_full_sim<- fixest::fepois(
    fml = fmla_full1,
    data = chopdat)
  # 4.2. Fit null model
  mod_null_sim <- fixest::fepois(
    fml = fmla_null1,
    data = chopdat)
  # 4.3. Extract the log-likelihood ratio test
  sim_result <- data.frame(
    rep = rep,
    lr_mod_ful_sim = logLik(mod_full_sim),
    lr_mod_null_sim = logLik(mod_null_sim),
    lr_test_null_sim = logLik(mod_full_sim) - logLik(mod_null_sim),
    seed = seed
  )
  
  return(sim_result)
}

## 2.2. Repeat the procedure lr_test_sim_null 10000 times ----------------------------
# to construct the null distribution of the test statistics

### Function to get independent seed for each run -------
get_seeds <- function(n_sim = 10000){
  set.seed(20251228, kind = "L'Ecuyer-CMRG", sample.kind = "Rejection")
  n_seed <- n_sim # Independent seed for each run in each scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}

### To repeat the procedure lr_test_sim_null 10000 times, we will use `foreach` loop,
#   which will be demonstrated in the example below (I have not found a way to make it into a function)

## 2.3. Calculate the p-value ------------------------------------------------------
#  for testing the null hypothesis that there does
#  not exist an interval with elevated risk for the outcome after vaccination

# Function to calculate p-value
mc_pvalue <- function(lrt_obs, lrt_null) {
  B <- nrow(lrt_null)
  (1 + sum(lrt_null$lr_test_sim >= lrt_obs)) / (1 + B)
}


# Demonstration: Applying the Xu 2013 method to itpdat datase ------------------

# This dataset is used in the paper Xu (2013)

## Step 1: Identify the optimal risk window

### Fit a series of fixed-effect models
risk_win <- seq(0, 364, by = 7)

itp_test_xu2013 <- foreach(k = seq_along(risk_win),
                           .combine = rbind) %do% 
  {
    i <- risk_win[k]
    
    data_analyse <- Xu_2013(
      formula = event ~ mmr + age,
      formula_null = event ~ age,
      indiv = case, 
      astart = sta,
      aend = end,
      aevent = itp,
      adrug = mmr,
      aedrug = mmr + i,
      expogrp = c(0),
      agegrp = c(427,488,549,610,671),
      data = itpdat,
      expo_of_interest = "mmr",
      risk_win_of_interest = 1
    )
    cbind(data_analyse, candidate_risk_win = i)
  }

### Identify the risk window corresponding to the highest LRT statistic
itp_maxlr <- xu2013_maxlr(result_table = itp_test_xu2013)
itp_maxlr

# The optimal risk window is 77 days, with the LRT statistics of 4.58.
# This is consistent with the results in Xu (2013) paper.

## Step 2: Obtain the p-value to test if there is indeed a risk window with elevated risk

### Obtain the distribution of the LRT statistic under the null
n_sim <- 10000
seeds <- get_seeds(n_sim) # Get 10,000 different seeds for 10,000 runs

# Set up parallel sessions to speed up the process
tic("Simulate null distribution")
plan(multisession, workers = 4)
lr_test_null_distribution <- foreach(i = 1:n_sim,
                                     .options.future = list(packages = c("SCCS", "dplyr", "survival"),
                                                            seed = TRUE),
                                     .combine = rbind # results of 10,000 run will be combined into a dataframe
) %dofuture% 
  { seed = seeds[i] # Assign the i_th seed for the i_th run
  one_lrt <- lr_test_sim_null(
    rep = i,
    seed = seed,
    maxlr = itp_maxlr,
    formula = event ~ vax_date_sim + age,
    formula_null = event ~ age,
    indiv = case, 
    astart = sta,
    aend = end,
    aevent = itp,
    adrug = vax_date_sim,
    aedrug = vax_date_sim + itp_maxlr$candidate_risk_win,
    expogrp = c(0),
    agegrp = c(427,488,549,610,671),
    data = itpdat,
    expo_of_interest = "mmr")
  
  one_lrt                                                            
  }
toc() # With parallel session on 4 cores, this procedure takes about 1 minute

### Calculate the p-value

p_val_itp <- mc_pvalue(
  lrt_obs  = itp_maxlr$lr_test,
  lrt_null = lr_test_null_distribution
)

# p-value = 0.0001, we reject the null hypothesis and conclude that there is indeed
# an elevated risk within the risk window of 77 days

# Example 2: dataset `condat` --------------------------------------------------
## Step 1: Identify the optimal risk window
### Fit a series of fixed-effect models
risk_win <- seq(7, 168, by = 7)
ageg <- seq(387,707,20)

con_test_xu2013 <- foreach(k = seq_along(risk_win),
                    .combine = rbind) %do% 
  {
    i <- risk_win[k]
    
    data_analyse <- Xu_2013(formula = event~hib+mmr+age,
                            formula_null = event ~ hib + age,
                            indiv = case, 
                            astart = sta,
                            aend = end,
                            aevent = conv,
                            adrug = cbind(hib, mmr),
                            aedrug = cbind(hib + 14, mmr + i),
                            expogrp = list(c(0), c(0)),
                            agegrp = ageg,
                            data = condat,
                            expo_of_interest = "mmr",
                            risk_win_of_interest = 1
    )
    cbind(data_analyse, candidate_risk_win = i)
    
  }

### Identify the risk window corresponding to the highest LRT statistic
con_maxlr <- xu2013_maxlr(result_table = con_test_xu2013)
con_maxlr

# The optimal risk window is 14 days, consistent with the one identified by the Xu 2011 method

## Step 2: Obtain the p-value to test if there is indeed a risk window with elevated risk

### Obtain the distribution of the LRT statistic under the null
n_sim <- 2000
seeds <- get_seeds(n_sim) 

# Set up parallel sessions to speed up the process
tic("Simulate null distribution")

plan(multisession, workers = 8)
con_lr_test_null_dist <- foreach(i = 1:n_sim, .options.future = list(packages = c("SCCS", "dplyr", "survival"),
                                                            seed = TRUE),
                      .combine = rbind # results of 10,000 run will be combined into a dataframe
) %dofuture% 
  { 
  seed = seeds[i] # Assign the i_th seed for the i_th run
  one_lrt <- lr_test_sim_null(
    rep = i,
    seed = seed,
    maxlr = con_maxlr,
    formula = event~hib+vax_date_sim+age,
    formula_null = event ~ hib + age,
    indiv = case, 
    astart = sta,
    aend = end,
    aevent = conv,
    adrug = cbind(hib, vax_date_sim),
    aedrug = cbind(hib + 14, vax_date_sim + con_maxlr$candidate_risk_win),
    expogrp = list(c(0), c(0)),
    agegrp = ageg,
    data = condat,
    expo_of_interest = "mmr")
  
  one_lrt                                                            
  }

toc() # 2000 runs take 5 mins with 8 cores

### Calculate the p-value

p_val_con <- mc_pvalue(
  lrt_obs  = con_maxlr$lr_test,
  lrt_null = con_lr_test_null_dist
)

#p-value = 0.0005
