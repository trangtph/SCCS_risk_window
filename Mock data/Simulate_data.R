set.seed(1246)

N <- 400

## ---- IDs and observation window ----
id_num  <- 175:(175 + N - 1)
obs_sta <- rep(1, N)
obs_end <- sample(400:577, N, replace = TRUE)

## ---- Exposure assignment ----
vax_group <- sample(c("bnt", "chad"), N, replace = TRUE)

bnt_1  <- rep(NA_real_, N)
bnt_2  <- rep(NA_real_, N)
chad_1 <- rep(NA_real_, N)
chad_2 <- rep(NA_real_, N)

for (i in 1:N) {
  if (vax_group[i] == "bnt") {
    bnt_1[i] <- sample(1:300, 1)
    if (runif(1) < 0.8 && bnt_1[i] + 30 <= obs_end[i]) {
      bnt_2[i] <- sample((bnt_1[i] + 30):obs_end[i], 1)
    }
  } else {
    chad_1[i] <- sample(1:300, 1)
    if (runif(1) < 0.8 && chad_1[i] + 30 <= obs_end[i]) {
      chad_2[i] <- sample((chad_1[i] + 30):obs_end[i], 1)
    }
  }
}

## ---- Event simulation ----
baseline_rate <- 1e-5
event_day <- rep(NA_real_, N)

for (i in 1:N) {
  
  days <- obs_sta[i]:obs_end[i]
  hazard <- rep(baseline_rate, length(days))
  
  # BNT risk windows
  for (t in c(bnt_1[i], bnt_2[i])) {
    if (!is.na(t)) {
      idx <- which(days >= t & days < t + 30)
      hazard[idx] <- hazard[idx] * 4
    }
  }
  
  # ChAd risk windows
  for (t in c(chad_1[i], chad_2[i])) {
    if (!is.na(t)) {
      idx <- which(days >= t & days < t + 30)
      hazard[idx] <- hazard[idx] * 2
    }
  }
  
  # simulate single event
  prob <- hazard / sum(hazard)
  event_day[i] <- sample(days, 1, prob = prob)
}

## ---- Assemble final dataset ----
sim_thrc_sccs <- data.frame(
  id = NA_character_,
  primary_diagnosis = NA_integer_,
  diagnosis_date = as.Date(NA),
  diagnosis_code = NA_character_,
  hosp_id_new = NA_character_,
  id_new = NA_character_,
  event_id = "THRC",
  event_name = NA_character_,
  obs_end_date = as.Date(NA),
  obs_sta_date = as.Date(NA),
  bnt_1_date = as.Date(NA),
  bnt_2_date = as.Date(NA),
  chad_1_date = as.Date(NA),
  chad_2_date = as.Date(NA),
  
  obs_sta = obs_sta,
  obs_end = obs_end,
  bnt_1  = bnt_1,
  bnt_2  = bnt_2,
  chad_1 = chad_1,
  chad_2 = chad_2,
  diagnosis = event_day,
  id_num = id_num
)

## ---- Write to CSV ----
write.csv(sim_thrc_sccs, here("Mock data", "simulated_thrc_sccs.csv"), row.names = FALSE)
