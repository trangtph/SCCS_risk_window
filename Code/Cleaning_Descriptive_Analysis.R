#############################################################
# Data Cleaning and Descriptive Analysis ---- ###############
#############################################################

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)
pacman::p_load(
  dplyr,
  rio,
  here,
  stringr,
  purrr,
  skimr
  )

# Load the tables of the CDM
comorbidity <- import(here("Mock data", "CDM", "comorbidity.sas7bdat"))
individual <- import(here("Mock data", "CDM", "individual.sas7bdat"))
outcome_event <- import(here("Mock data", "CDM", "outcome_event.sas7bdat"))
covid_infection <- import(here("Mock data", "CDM", "covid_infection2.csv"))
covid_vacc <- import(here("Mock data", "CDM", "covid_vacc2.csv"))


# Check the structure and format of each table ----
cdm <- list(individual = individual, 
            covid_infection = covid_infection, 
            covid_vacc = covid_vacc,
            outcome_event = outcome_event, 
            comorbidity = comorbidity)

lapply(cdm, skimr::skim)

## Format variables ----
### Clean character missingness
clean_missing <- function(df) {
  char_cols <- sapply(df, is.character)
  
  df[char_cols] <- lapply(df[char_cols], function(x) {
    x <- trimws(x) # Remove leading and/or trailing whitespace from character strings
    x[x %in% c("", "NA")] <- NA
    x
  })
  df
}

### Convert date columns to date format
convert_date <- function(df, date_pattern = "_date(_\\d+)?$", #"_date" optionally followed by "_" and one or more digits
                        format = "%d-%m-%Y") {
  
  date_cols <- grep(date_pattern, names(df), value = TRUE) #get the names of all the columns containing "_date"
  
  for (col in date_cols) {
    x <- df[[col]]
    
    if (!is.character(x)) next
    
    parsed <- as.Date(x, format = format)
    
    # Warn if conversion failed for non-missing values
    failed <- !is.na(x) & is.na(parsed)
    if (any(failed)) {
      warning(
        sprintf(
          "Date parsing failed in column '%s' for %d rows",
          col, sum(failed)
        ),
        call. = FALSE
      )
    }
    df[[col]] <- parsed
  }
  df
}

cdm_cleaned <- lapply(cdm, function(df){
  df |> clean_missing() |>
    convert_date()
})

## Add calculated variables

cdm$individual$age <- 2020 - cdm$individual$birth_year

# Adding names of the outcome and covariate ----

 <- read_excel(here("", 
                                 ".xlsx"), 
                            sheet ="")

# 1. Apply exclusion criteria -----------------------------------------------------

## General exclusion criteria ------------------------
### Age 12 years 

id_12 <- cdm$individual %>% 
  dplyr::filter(age < 12) %>% 
  dplyr::distinct(id)

### obs period 1st Oct 2020 - 30 Apr 2022
outcome_before_obs <-  cdm$outcome_event %>% 
  dplyr::filter(!is.na(diagnosis_date) & diagnosis_date < as.Date("01-10-2020", format = "%d-%m-%Y")) %>% 
  dplyr::distinct(id, event_name)

outcome_after_end <- cdm$outcome_event %>% 
  dplyr::filter(!is.na(diagnosis_date) & diagnosis_date > as.Date("30-04-2022", format = "%d-%m-%Y")) %>% 
  dplyr::distinct(id, event_name)

### Exclude missing year of birth, outcome date, COVID-19 vacc date, vaccine brand or dose



### Exclude heterologious dose (For BNT & ChAd only)

### Exclude individuals with <3 years of data availability prior to 1st Oct 2020

## Specific exclusion criteria to each outcome -------
### Create three tables for three outcomes

### Exclude Inpatient diagnosis of the outcome between 1st October 2017 and 30th September 2020

### Only keep the first event (for each outcome) per individual

## Final list of ids for the whole sample and per outcome --------

# 2. Extract censoring information ---------------------------------------------

## 1st dose of COVID-19 vaccine other than BNT-162b2 or ChAdOx1

## 3rd dose of BNT-162b2 or ChAdOx1

## 1st COVID-19 infection

## Merge information to table `individual` and pick the earliest date before 30 Apr 2022 as censoring date

# 3. Descriptive analysis ------------------------------------------------------

## Number & distribution of outcome events --------

## Distribution of vacc date ----

## Distribution of 1st infection ----

## Check assumptions of SCCS model ----

### Event-dependent observation period: histogram of time from event until end of observation

### Event-dependent exposure: time interval between the start of each exposure and the event

## Baseline characteristics of cases (per outcome) ----

# 4. Build exposure-outcome tables in SCCS-compatible format -------------------










