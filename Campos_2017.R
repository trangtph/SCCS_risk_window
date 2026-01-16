#############################################################
# Analysis code for Campos (2017) ---------------------------
#############################################################

# 1. Setup Environment and Load Data --------------------------
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
  scales,
  splines2,
  rlang 
)

# Define functions ----------
# First part on SCCS function definition is similar to the script for Xu method (2011) - by Trang

# Function Campos_2017 performs the following steps:
# (1) Data preprocessing and reshaping
#     - Processes input variables and reformats the data into SCCS-compatible
#       interval format using SCCS::formatdata.
#     - This follows the same internal procedure as in `SCCS::standardsccs`.
# (2) SCCS model fitting and inference
#     - Fits a conditional Poisson (clogit) SCCS model with individual-level
#       stratification and interval offsets.
#     - Extracts effect estimates corresponding to the specified risk window,
#       including IRR, standard error, p-value, and 95% confidence interval.
# Function Campos_plot performs the following:
# (3) takes the result from Campos_2017 function and makes a scatterplot/line function 
#     using knot candidates
# (4) calculate SSE and find the know with minimum SSE. Knot with minimum SSE corresponds 
#     to the optimal day length
# (5) plot to visualise


# 2. Define Function Campos_2017 -------------------------------

# Function Campos_2017 performs data preprocessing, reshaping, and SCCS model fitting
Campos_2017 <- function(formula, #The dependent variable should always be "event"
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
  
  
  # Fit SCCS Model
  fmla <- paste(formula, "+", "strata(indivL)", "+", "offset(log(interval))")
  fmla1 <- as.formula(paste("event~", fmla[3]))
  mod <- clogit(formula = fmla1, data = chopdat)
  sum_mod <- summary(mod)
  # Export model statistics
  main_expo <- paste0(expo_of_interest, risk_win_of_interest)
  est_L <- stats::coef(sum_mod)[main_expo, 'coef']
  IRR_L <- stats::coef(sum_mod)[main_expo, 'exp(coef)']
  se_L <- stats::coef(sum_mod)[main_expo, 'se(coef)']
  p_val <- stats::coef(sum_mod)[main_expo, 'Pr(>|z|)']
  IRR_L_low_CI <- sum_mod$conf.int[main_expo, 'lower .95']
  IRR_L_up_CI <- sum_mod$conf.int[main_expo, 'upper .95']
  
  model_stat <- data.frame(est_L, IRR_L,
                           se_L, p_val,
                           IRR_L_low_CI, IRR_L_up_CI)
  
  return(model_stat)
  
}

# 3. Define the Campos_plot Function -----------------

# Function `Campos_plot` takes the results from `Campos_2017` function, 
# and makes a scatterplot/line plot using integer days (e.g. 21-147) as potential knots.

##- adjusted = define knot in DAYS

Campos_plot <- function(data, 
                        risk_win,# variable 'risk length' in the data
                        T_L_inv, # variable '1/risk length' in the data
                        IRR_L,   # variable 'estimated IRR' in the data
                        xbreak = 20 # number of breaks for the x axis plot
                        ){ 
  
  # 1. Define range of candidate knots
  #day_candidates <- seq(from = start_day, to = end_day, by = 1)
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
  
  # 4. Find the corresponding optimal day length (t = 1 / (1/t))
  candidate_risk_win <- data[[risk_win]]
  best_day_length <- candidate_risk_win[which.min(sse_values)]
  
  # 5. Fit final spline using the best knot
  lp  <- pmax(0, x - best_knot_inv)
  lp2 <- lp^2
  
  final_model <- lm(y ~ x + lp + lp2)
  
  # 6. Generate prediction grid for plotting
  x_grid <- seq(min(x), max(x), length.out = 200)
  
  grid_df <- data.frame(
    x   = x_grid,
    lp  = pmax(0, x_grid - best_knot_inv),
    lp2 = pmax(0, x_grid - best_knot_inv)^2
  )
  
  grid_df$pred_y <- predict(final_model, newdata = grid_df)
  
  # 7. Visualise the results
  optimal_point <- data[data[[T_L_inv]] == best_knot_inv, ]
  
  p <- ggplot(data, aes(x = .data[[T_L_inv]], y = .data[[IRR_L]])) +
    geom_point(alpha = 0.8, color = "black") +
    geom_line(data = grid_df, aes(x = x, y = pred_y), linewidth = 1) +
    geom_point(data = optimal_point, 
               aes(x = .data[[T_L_inv]], .data[[IRR_L]]),
               color = "darkred", 
               size = 3) + 
    labs(title = "Optimal Quadratic Spline Fit (Campos 2017)",
         subtitle = paste("Optimal knot at", best_day_length, "days"), 
         y = "IRR (L)") +
    scale_x_continuous(
      name = "Inverse Risk Length (1/T(L))",
      breaks = scales::breaks_pretty(n = xbreak)
    ) +
    theme_minimal()
  
  return(list(model = final_model, best_knot_inv = best_knot_inv, 
              best_day_length = best_day_length, plot = p))
}


# 4. Execute the Analysis --------------------------------------

# Example 1: Use 'itpdat' 
# Load the example data from the SCCS package
data(itpdat)


# Define the risk windows used for generating the scatterplot points (7-day increments)
risk_win <- seq(21, 147, by = 7)  # follow the same format as described in the paper

itp_test_campos <- foreach(k = seq_along(risk_win),
                           .combine = rbind) %do% 
  {
    i <- risk_win[k]
    
    res <- Campos_2017(event ~ mmr + age,
                       indiv = case, 
                       astart = sta,
                       aend = end,
                       aevent = itp,
                       adrug = mmr,
                       aedrug = mmr + i,
                       expogrp = c(0,1),
                       agegrp = c(427, 488, 549, 610, 671),
                       data = itpdat,
                       expo_of_interest = "mmr",
                       risk_win_of_interest = 2)
    
    # Add required columns for plotting: risk length (t) and inverse (1/t)
    res$candidate_risk_win = i
    res$inverse_t <- 1/i
    res
  }

# View the collected results data frame
View(itp_test_campos)

# Use the plotting function to find the optimal knot from integer day candidates
results <- Campos_plot(
  data = itp_test_campos,
  risk_win = "candidate_risk_win",
  T_L_inv = "inverse_t",
  IRR_L = "IRR_L"
)

# View the generated plot
print(results$plot)
# Print the IRR that corresponds to the optimal risk length
print(itp_test_campos[itp_test_campos$candidate_risk_win==results$best_day_length,])


## Results are consistent with Campos (2017 paper), tau = 84 days