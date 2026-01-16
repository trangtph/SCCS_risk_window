#############################################################
# Analysis code for Xu (2011) --------------- ###############
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
  tictoc
)

options(scipen = 999)

# 1. Define functions -------------------------------------------------------------

## 1.1. Function to fit a series of SCCS models ------------------------------------
# Function Xu_2011 performs the following steps:
# (1) Data preprocessing and reshaping
#     - Processes input variables and reformats the data into SCCS-compatible
#       interval format using SCCS::formatdata.
#     - This follows the same internal procedure as in `SCCS::standardsccs`.
# (2) Average time at risk calculation
#     - For a given risk window L, computes the average time at risk:
#           T(L) = (1 / n) * sum_{i,j} t_{ij1}
#       where:
#         i indexes individuals,
#         j indexes age (or season) intervals,
#         t_{ij1} is the time spent in the risk window of interest,
#         n is the number of individuals contributing to that risk window.
# (3) SCCS model fitting and inference
#     - Fits a conditional Poisson (clogit) SCCS model with individual-level
#       stratification and interval offsets.
#     - Extracts effect estimates corresponding to the specified risk window,
#       including IRR, standard error, p-value, and 95% confidence interval.

# The way the arguments in the function Xu_2011 are specified is exactly like how
# the arguments in `SCCS::standardsccs` are specified.
# Two additional arguments are `expo_of_interest` and `risk_win_of_interest`.
#   `expo_of_interest`: Character string giving the name of the exposure variable 
#                       for which the risk window is being scanned.
#   `risk_win_of_interest`: Numeric index indicating which exposure level 
#                           (as defined by `expogrp`) corresponds to the risk window of interest.
#
#     Example:
#       - If expogrp = 0:
#           There are two exposure levels: 0 = control, 1 = risk window (day 0 to L)
#           To scan over window lengths L, set risk_win_of_interest = 1.
#       - If expogrp = c(0, 1):
#           There are three exposure levels:0 = control, 1 = day 0 only, 2 = days 1 to L
#           To scan over window lengths L, set risk_win_of_interest = 2.


Xu_2011 <- function(formula, #The dependent variable should always be "event"
                    indiv, astart, aend, aevent, adrug, aedrug, expogrp = list(), washout = list(), 
                    sameexpopar = list(), agegrp = NULL, seasongrp=NULL, dob=NULL, dataformat="stack", data,
                    expo_of_interest, risk_win_of_interest = 2 # 0 = control, 1 = day 0
                          ) 
  {
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
  
  
  # ------------------Getting the fixed covariates from the formula ------------#
  
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
  
  
  expo_dat <- chopdat[chopdat[[expo_of_interest]] == risk_win_of_interest, ]
  T_L <- sum(expo_dat$interval) / length(unique(expo_dat$indivL))
  T_L1 <- 1/ T_L
  
  # Fit SCCS Model
  fmla <- paste(formula, "+", "strata(indivL)", "+", "offset(log(interval))")
  fmla1 <- as.formula(paste("event~", fmla[3]))
  mod <- clogit(formula = fmla1, data = chopdat)
  sum_mod <- summary(mod)
  # Export model statistics
  main_expo <- paste0(expo_of_interest, risk_win_of_interest)
  est_L <- stats::coef(sum_mod)[main_expo, 'coef'] # 2nd row because 1st row is for day 0.
  IRR_L <- stats::coef(sum_mod)[main_expo, 'exp(coef)']
  se_L <- stats::coef(sum_mod)[main_expo, 'se(coef)']
  p_val <- stats::coef(sum_mod)[main_expo, 'Pr(>|z|)']
  IRR_L_low_CI <- sum_mod$conf.int[main_expo, 'lower .95']
  IRR_L_up_CI <- sum_mod$conf.int[main_expo, 'upper .95']
  
  model_stat <- data.frame(T_L, T_L1, est_L, IRR_L,
                           se_L, p_val,
                           IRR_L_low_CI, IRR_L_up_CI)
  
  return(model_stat)
  
}

## 1.2. Function to visualise the relationship between IRR_L and 1/T(L) --------
# Function `Xu_2011_plot` takes the results from `Xu_2011` function, and make a 
# scatterplot and line plot of IRR_L versus 1/T(L) (average time at risk).
# The maximum IRR is highlighed on the plot.

# To interpret the plot: the change-point marking the start of the linear region 
# is the optimal risk window. Most of the time it also yields the maximum IRR 
# (not true if event is sparse).

Xu_2011_plot <- function(data,
                         x = "T_L1",
                         y = "IRR_L",
                         xlab = "1/ T(L)",
                         ylab = "R (L)", 
                         xbreak = 20) {
  
  
  data2 <- data %>% 
    filter(se_L < 10)
  
  # row with maximum IRR
  max_row <- data2 %>% 
    slice_max(.data[[y]], n = 1, with_ties = FALSE)
  
plot_irr <- ggplot(data2, aes(x = .data[[x]], y = .data[[y]])) +
    geom_line() +
    geom_point() +
    geom_point(
      data = max_row,
      size = 3,
      colour = "darkred"
    ) +
    scale_x_continuous(
      name = xlab,
      breaks = scales::breaks_pretty(n = xbreak)
    ) +
    labs(y = ylab) +
  annotate(
    "text",
    x = max_row[[x]],
    y = max_row[[y]],
    label = paste0("Max IRR at L = ", round(max_row$candidate_risk_win, 2)),
    vjust = -1,
    hjust = 0.5
  ) +
    theme_minimal()

plot_irr
}

# 2. Test the function on datasets from the `SCCS` package ---------------------------------------------------

## 2.1. Test for dataset 'condat' --------------------------------------------------

# Note: this dataset is similar to our COVID-19 vaccine dataset with multiple exposures.
# When scan the risk window for e.g. mmr, we fix the risk window of the other exposure
# I specify `expogrp = list(c(0,1), c(0,1)),` because we want to isolate day 0, and are
# interested in scanning the risk window starting on day 1.

risk_win <- seq(7, 168, by = 7)
ageg <- seq(387,707,20)

con_test <- foreach(k = seq_along(risk_win),
                    .combine = rbind) %do% 
  {
    
    i <- risk_win[k]
    
    data_analyse <- Xu_2011(event~hib+mmr+age,
                            indiv = case, 
                            astart = sta,
                            aend = end,
                            aevent = conv,
                            adrug = cbind(hib, mmr),
                            aedrug = cbind(hib + 14, mmr + i),
                            expogrp = list(c(1), c(1)),
                            agegrp = ageg,
                            data = condat,
                            expo_of_interest = "mmr",
                            risk_win_of_interest = 1
    )
    cbind(data_analyse, candidate_risk_win = i)
    
  }

Xu_2011_plot(data = con_test)


## 2.2. Test for dataset 'itpdat' --------------------------------------------------

# This dataset was used for as example in the paper Xu et al. (2011)

risk_win <- seq(21, 147, by = 7)

itp_test <- foreach(k = seq_along(risk_win),
                    .combine = rbind) %do% 
  {
    
    i <- risk_win[k]
    
    data_analyse <- Xu_2011(event~mmr+age,
                            indiv = case, 
                            astart = sta,
                            aend = end,
                            aevent = itp,
                            adrug = mmr,
                            aedrug = mmr + i,
                            expogrp = c(0,1),
                            agegrp = c(427,488,549,610,671),
                            data = itpdat,
                            expo_of_interest = "mmr",
                            risk_win_of_interest = 2
    )
    cbind(data_analyse, candidate_risk_win = i)
    
  }

Xu_2011_plot(data = itp_test)

# Comment: Max IRR at L_M = 35, but there is no linear relationship when L > L_M (1/T(L) < 1/T(L_M)). 
# We choose the 2nd highest IRR of 3.3 at L_M = 77 day because it shows linear relationship when L > 77. 
# This is consistent with the results of this example in the paper Xu et al. (2011)

## 2.3. Test for dataset 'gbsdat' -------------------------------------------------

seas <- cumsum(c(31,30,31,31,28,31)) # calendar month

risk_win <- seq(1, 70, by = 1)

gbs_test2 <- foreach(k = seq_along(risk_win),
                     .combine = rbind) %do% {
                       
                       i <- risk_win[k]
                       
                       data_analyse <- Xu_2011(event~flu+age,
                                               indiv = case, 
                                               astart = sta,
                                               aend = end,
                                               aevent = gbs,
                                               adrug = flu,
                                               aedrug = flu + i,
                                               expogrp = c(0,1),
                                               agegrp = seas,
                                               data = gbsdat,
                                               expo_of_interest = "flu",
                                               risk_win_of_interest = 2
                       )
                       cbind(data_analyse, candidate_risk_win = i)
                       
                     }

Xu_2011_plot(data = gbs_test2)

