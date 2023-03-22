# Comparison of functionality for datadelay -------------------------------


# Load libraries and data -------------------------------------------------

# Install libraries
library(tidyverse)
library(incidence2)

#setwd("~/Documents/GitHub/epiverse-trace")
install("datadelay")
library(datadelay)

remotes::install_github("epiverse-trace/epiparameter")
library(epiparameter)

remotes::install_github("epiforecasts/EpiNow2@develop")
library(EpiNow2)

library(coarseDataTools)

# Define onset-to-death
onset_to_death_ebola <- epiparameter::epidist_db(
  disease = "Ebola Virus Disease", 
  epi_dist = "onset_to_death",
  author = "Barry_etal")

# Set ebola-cfr directory
#setwd("~/Documents/GitHub/ebola-cfr/")

# Load line list 
data_ebola <- read_csv("data/ebola_1976_linelist.csv")

cut_off <- as.Date("1976-10-05") # date to analyse up to in real-time

# Calculate CFR based on line list and ground truth ------------------------

df_ll_subset_onset <- subset(data_ebola, disease_onset <= cut_off) # onsets up to cut off

df_ll_subset_outcome <- subset(data_ebola, disease_ended <= cut_off) # outcomes up to cut off

# Estimate CFR based on known outcomes in real-time
CFR_binomial_RT_ll <- prop.test(x=sum(df_ll_subset_outcome$status=="died"), n=nrow(df_ll_subset_outcome))
dt_CFR_RT_ll <- c(CFR_binomial_RT$estimate,CFR_binomial_RT$conf.int) |> format_cfr_neatly(type = "Line list CFR (RT)")

# Estimate CFR based on future follow ups of known onsets
CFR_binomial_ll_true <- prop.test(x=sum(df_ll_subset_onset$status=="died"), n=nrow(df_ll_subset_onset))
dt_CFR_true_ll <- c(CFR_binomial_ll_true$estimate,CFR_binomial_ll_true$conf.int) |> format_cfr_neatly(type = "Line list CFR (true hindsight)")


# Convert Ebola line list into incidence ----------------------------------
# Convert cases and deaths to incidence series
case_data <- incidence2::incidence(data_ebola,disease_onset) |> complete_counts() |> rename(cases = count,date=date_index)

data_ebola_died <- data_ebola |> filter(status=="died")
death_data <- incidence2::incidence(data_ebola_died,disease_ended) |> complete_counts() |> rename(deaths = count,date=date_index)

data_ebola_ts <- merge(case_data,death_data,all=T) # merge timeseries
data_ebola_ts[is.na(data_ebola_ts)] <- 0 # replace NA with zero counts

# Define subset of data for analysis
df_ebola_subset <- subset(data_ebola_ts, date <= cut_off)

# Plot data
par(mfcol=c(1,1),mar=c(3,4,1,1),las=1)
plot(df_ebola_subset$date,df_ebola_subset$cases,
     type="s",lwd=2,col="black",ylab="cases (black) and deaths (orange)",xlab="",ylim=c(0,20),yaxs="i",main="Ebola 1976") 
lines(df_ebola_subset$date,df_ebola_subset$deaths,type="s",lwd=2,col="orange")



# Estimate CFR with datadelay -------------------------------------------------
dt_ncfr_static_ebola <- static_cfr(
  df_ebola_subset,
  correct_for_delays = FALSE) |> 
  format_cfr_neatly(type = "Naive")

dt_ccfr_static_ebola <- static_cfr(
  df_ebola_subset,
  correct_for_delays = TRUE,
  epi_dist = onset_to_death_ebola) |>
  format_cfr_neatly(type = "Corrected")


# Estimate with EpiNow2 ---------------------------------------------------

# Format parameters
onset_to_death_meansd <- epiparameter::gamma_shapescale2meansd(
  shape = 2.4, # NOTE HARD CODING - NEED TO UPDATE WITH NEW EPIPARAMETER FUNCTIONALITY WHEN READY
  scale = 1/0.3 # NOTE HARD CODING
)
onset_to_death_logmean <- EpiNow2::convert_to_logmean(
  mean = onset_to_death_meansd[["mean"]],
  sd = onset_to_death_meansd[["sd"]]
)
onset_to_death_logsd <- EpiNow2::convert_to_logsd(
  mean = onset_to_death_meansd[["mean"]], # NEED TO CHECK CONVERSION FROM GAMMA
  sd = onset_to_death_meansd[["sd"]]
)

# Format data
data_epinow2 <- df_ebola_subset |> rename(primary = cases, secondary = deaths)

cases_to_deaths <- EpiNow2::estimate_secondary(
  data_epinow2,
  delays = delay_opts(list(
    mean = onset_to_death_logmean,
    mean_sd = 0.1,
    sd = onset_to_death_logsd,
    sd_sd = 0.1,
    max = 21)
  ),
  secondary = secondary_opts(type = "incidence"),
  obs = obs_opts(
    scale = list(mean = 0.5, sd = 10),
    family = "poisson",
    week_effect = FALSE
  ),
  verbose = FALSE
)

# Plot outputs
cfr_samples <- rstan::extract(cases_to_deaths$fit, "frac_obs")[[1]][, 1]

tibble(cfr = cfr_samples) |>
  ggplot(aes(x = cfr)) +
  geom_density() +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.25) +
  theme_bw() +
  xlab("CFR") +
  ylab("Posterior density")

dt_ccfr_static_ebola_epinow2 <- quantile(cfr_samples, c(0.5, 0.025, 0.975)) |> format_cfr_neatly(type = "EpiNow2 Corrected") ## median + 95% CI

# Compare datadelay and EpiNow2
dt_static_clean_ebola <- rbind(dt_ncfr_static_ebola, # Naive incidence calculation
                               dt_ccfr_static_ebola, # datadelay adjusted estimate
                               dt_ccfr_static_ebola_epinow2, # EpiNow2 estimate
                               dt_CFR_RT_ll, # Based on known outcomes in real-time
                               dt_CFR_true_ll # Based on subsequent known outcomes of real-time onsets (i.e. ground truth)
                               )

