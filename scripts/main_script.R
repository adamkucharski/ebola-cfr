# - - - - - - - - - - - - - - - - - - - - - - - 
# CFR estimation analysis
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load libraries
library(tidyverse)
library(devtools)
library(incidence2)
library(coarseDataTools)
install_github("epiverse-trace/epiparameter")
library(epiparameter)
library("EpiNow2") ## remotes::install_github("epiforecasts/EpiNow2")
library("rstan")

# Define probability mass function for onset-to-death
onset_to_death_ebola <- epiparameter::epidist("ebola","onset_to_death")$pmf
onset_to_death_param <- epiparameter::epidist("ebola","onset_to_death")$param
onset_to_death_meansd <- epiparameter::gamma_shapescale2meansd(
  shape = onset_to_death_param[["shape"]],
  scale = onset_to_death_param[["scale"]]
)
onset_to_death_logmean <- EpiNow2::convert_to_logmean(
  mean = onset_to_death_meansd[["mean"]],
  sd = onset_to_death_meansd[["sd"]]
)
onset_to_death_logsd <- EpiNow2::convert_to_logsd(
  mean = onset_to_death_meansd[["mean"]],
  sd = onset_to_death_meansd[["sd"]]
)

xx <- 0:20
plot(xx,onset_to_death_ebola(xx),xlab="days")

# Load data and functions
# setwd("~/Documents/GitHub/ebola-cfr/")
source("R/data_load.R")
source("R/cfr_function.R")

# Run example delay adjusted CFR estimation on 1976 Ebola data: -----------
cutoff <- 30
data_in <- data1976 |> # Incidence data
  slice(seq_len(cutoff))
data_epinow2 <- data_in |>
  rename(primary = cases, secondary = deaths)

correct_cfr <- correct_cfr_1976 # Actual final CFR

# Estimate CFR over first 30 days
out_cfr <- scale_cfr(case_incidence = data_in$cases,
                     death_incidence = data_in$deaths,
                     delay_fun = onset_to_death_ebola)

# Unadjusted ('naive') CFR and 95% CI:
out_cfr$nCFR

# Delay adjusted CFR and 95% CI:
out_cfr$cCFR

# EMforCFR method in coarseDataTools (Note: work in progress):
## reporting.param <- rep(0,nrow(data1976_all_pick)-1)
## out_em_cfr <- EMforCFR(assumed.nu = onset_to_death_ebola(xx), alpha.start.values = reporting.param, full.data = data1976_all_pick)
## out_em_cfr <- EMforCFR(assumed.nu = onset_to_death_ebola(xx), alpha.start.values = reporting.param, full.data = data_cdt)


# Plot estimates over time on both 1976 and 2014 data ------------------------------------------

# Black lines, daily new cases; orange lines, daily new deaths
# Red lines, naive CFR calculation, i.e. deaths/cases, with 95% CI; blue lines, delay-adjusted calculation with 95%

source("R/plot_examples.R")

cases_to_deaths <- estimate_secondary(
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

cfr_samples <- rstan::extract(cases_to_deaths$fit, "frac_obs")[[1]][, 1]
quantile(cfr_samples, c(0.025, 0.5, 0.975)) ## median + 95% CI
tibble(cfr = cfr_samples) |>
  ggplot(aes(x = cfr)) +
  geom_density() +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.25) +
  theme_bw() +
  xlab("CFR") +
  ylab("Posterior density")


## forecast from the following 30 days of case data (this could be replaced by a case forecast)
cases <- data1976 |>
  slice(seq(cutoff + 1, cutoff * 2)) |>
  select(date, value = cases)
forecast <- forecast_secondary(cases_to_deaths, cases)
plot(forecast)
