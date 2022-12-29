# - - - - - - - - - - - - - - - - - - - - - - - 
# Comparison of EpiLine and Dynamical Truncation
# Examples from below repos, reproduced here for illustrative comparison:
# https://github.com/parksw3/dynamicaltruncation/
# https://github.com/BDI-pathogens/EpiLine
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load libraries
library(tidyverse)
library(devtools)
library(data.table) # data.table::update_dev_pkg()
library(purrr)
#library(rstan)
#library(brms)

install_github("BDI-pathogens/EpiLine")
library(EpiLine)

install_github("parksw3/dynamicaltruncation")
library(dynamicaltruncation)

# EpiLine simulation ------------------------------------------------------

set.seed( 1 )

# define the length of the simulatiopn
t_rep          <- 50 # length of time for which data is reported
t_symptom_pre  <- 30 # time before the reporting period to simulate
t_symptom_post <- 5  # time after the reporting period to simulate
t_max          <- t_rep + t_symptom_post + t_symptom_pre

# set up the variable r(t) and distribution
symptom_0 <- 2                                # initial number of symptomatic people
r         <- 0.1 - 0.13 * ( 1:t_max ) / t_max # r(t) in the simulation
xi        <- -1 + 6 * ( t_max:1 ) / t_max          # xi parameter in the symptom-report dist
lambda    <- 2 + ( t_max:1 ) / t_max         # lambda parameter in the symptom-report dist

simulation <- symptom_report.simulator(
  t_rep          = t_rep,
  t_symptom_pre  = t_symptom_pre,
  t_symptom_post = t_symptom_post,
  symptom_0   = symptom_0,
  r           = r,
  dist_xi     = xi,
  dist_lambda = lambda
)

# Fit model:

# data 
reported    <- simulation$reported
ll_report   <- simulation$linelist$report
ll_symptom  <- simulation$linelist$symptom
report_date <- as.Date("2022-04-01")

# fit using model
mcmc_n_samples <- 100
mcmc_n_chains  <- 1
fit <- symptom_report.fit( reported, ll_symptom, ll_report, report_date = report_date, 
                           mcmc_n_samples = mcmc_n_samples, mcmc_n_chains = mcmc_n_chains )

fit$plot.symptoms( simulation = simulation)

# Estimation with dynamical truncation ------------------------------------

outbreak <- simulate_gillespie(seed = 101)
secondary_dist <- data.table( meanlog = 1.8, sdlog = 0.5) |> add_natural_scale_mean_sd()

obs <- outbreak |>
  simulate_secondary(
    meanlog = secondary_dist$meanlog[[1]],
    sdlog = secondary_dist$sdlog[[1]]
  ) |> observe_process()

truncated_obs <- obs |>
  filter_obs_by_obs_time(obs_time = 25) |>
  DT(sample(1:.N, 200, replace = FALSE))

truncated_cases <- construct_cases_by_obs_window(
  obs, windows = c(25), obs_type = "stime"
)

plot_cases_by_obs_window(truncated_cases)

naive_fit <- naive_delay(data = truncated_obs, cores = 1, refresh = 0)


