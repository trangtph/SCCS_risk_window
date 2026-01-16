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

Campos_plot <- function(data, x_var_inv, y_var, start_day, end_day) {
  
  # 1. Define range of possible knots in DAYS (t = 21, ..., 147)
  # We use +1 and -1 to ensure candidates are inside the data range [1/147, 1/21]
  day_candidates <- seq(from = start_day + 1, to = end_day - 1, by = 1)
  # Convert to the inverse scale (1/t) for the x-axis variable
  knot_candidates_inv <- 1 / day_candidates
  
  # Get data boundaries for fixed splines
  x_data <- data[[x_var_inv]]
  b_knots <- range(x_data, na.rm = TRUE)
  
  # 2. Define internal SSE calculation function (operates on INVERSE values)
  calc_sse <- function(k_inv, df_internal, x_inv, y, boundaries) {
    # Ensure knots are strictly inside boundaries to avoid the error
    fit <- lm(get(y) ~ splines2::bSpline(get(x_inv), knots = k_inv, 
                                         Boundary.knots = boundaries, degree = 2), 
              data = df_internal)
    sum(resid(fit)^2)
  }
  
  # 3. Find the knot (in inverse scale) that minimizes SSE
  sse_values <- sapply(knot_candidates_inv, calc_sse, df_internal = data, 
                       x_inv = x_var_inv, y = y_var, boundaries = b_knots)
  best_knot_inv <- knot_candidates_inv[which.min(sse_values)]
  
  # 4. Find the corresponding optimal day length (t = 1 / (1/t))
  best_day_length <- 1 / best_knot_inv
  
  # 5. Fit final model using the best inverse knot location
  final_fmla <- as.formula(paste0(y_var, " ~ splines2::bSpline(", x_var_inv, ", knots = ", best_knot_inv, ", degree = 2)"))
  final_model <- lm(final_fmla, data = data)
  
  # 6. Generate prediction grid for smooth plotting
  t_grid <- data.frame(seq(b_knots[1], b_knots[2], length.out = 200))
  colnames(t_grid) <- x_var_inv
  t_grid$pred_y <- predict(final_model, newdata = t_grid)
  
  # 7. Create the plot using .data pronoun for updated ggplot compatibility
  p <- ggplot(data, aes(x = .data[[x_var_inv]], y = .data[[y_var]])) +
    geom_point(alpha = 0.6, color = "darkgray") +
    geom_line(data = t_grid, aes(x = .data[[x_var_inv]], y = .data[["pred_y"]]), 
              color = "blue", linewidth = 1) +
    geom_vline(xintercept = best_knot_inv, linetype = "dashed", color = "red") +
    labs(title = "Optimal Quadratic Spline Fit (Campos 2017)",
         subtitle = paste("Optimal knot at", round(best_day_length, 1), "days"),
         x = "Inverse Risk Length (1/t)", y = "Relative Risk") +
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
    res$risk_length <- i
    res$inv_t <- 1/i
    res
  }

# View the collected results data frame
View(itp_test_campos)

# Use the plotting function to find the optimal knot from integer day candidates
results <- Campos_plot(
  data = itp_test_campos, 
  x_var_inv = "inv_t", 
  y_var = "IRR_L",
  start_day = 21,
  end_day = 147
)

# Print the final optimal day length and view the generated ggplot object
cat("The optimal risk length calculated from integer day candidates is:", results$best_day_length, "days.\n")
print(results$plot)
# print the IRR that corresponds to the optimal risk length
print(itp_test_campos$IRR_L[itp_test_campos$risk_length==results$best_day_length])


## Note: optimal knot at 0.01301546 >> corresponding to tau = 77 days
## not consistent with results described in paper, tau = 84 days