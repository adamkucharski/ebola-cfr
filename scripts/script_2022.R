# - - - - - - - - - - - - - - - - - - - - - - - 
# CFR estimation analysis
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load libraries
library(tidyverse)
library(devtools)
library(incidence2)
install_github("epiverse-trace/epiparameter")
library(epiparameter)
library(fitdistrplus)

# setwd("~/Documents/GitHub/ebola-cfr/")

# Define probability mass function for onset-to-death
onset_to_death_ebola <- epiparameter::epidist("ebola","onset_to_death")
onset_to_death_ebola <- onset_to_death_ebola$pmf # get pmf

xx <- 0:20
plot(xx,onset_to_death_ebola(xx),xlab="days")

# Load data and functions
source("R/cfr_function.R")


# Import 2022 Ebola data --------------------------------------------------
data_ebola <- read_csv("https://3mmuwilir3.execute-api.eu-central-1.amazonaws.com/web/url?folder=&file_name=latest.csv")
# https://raw.githubusercontent.com/globaldothealth/ebola/main/latest.csv

# Estimate CFR for cases with known outcomes --------------------------------------------------
data_ebola <- data_ebola |> mutate(Date_case = as.Date(ifelse(!is.na(Date_onset),Date_onset,Date_confirmation),"1970-01-01"))
data_outcome <- data_ebola %>% filter(!is.na(Date_Death) | !is.na(Date_Recovered))


# Estimate CFR from incidence data --------------------------------------------------
case_data <- incidence2::incidence(data_ebola,Date_case) |> complete_counts()
death_data <- incidence2::incidence(data_ebola,Date_Death) |> complete_counts()

# Plot data:
plot(case_data$date_index,case_data$count,ylim=c(0,1.2*max(case_data$count)),
     type="s",lwd=2,col="black",ylab="cases and deaths",xlab="",yaxs="i",main="") # Plot naive CFR
lines(death_data$date_index,death_data$count,type="s",lwd=2,col="orange") # Plot naive CFR

# Estimate CFR:
out_cfr <- scale_cfr(case_incidence = case_data$count,
                     death_incidence = death_data$count,
                     delay_fun = onset_to_death_ebola)


# Output estimate based on known outcomes
CFR_binomial <- prop.test(x=sum(!is.na(data_outcome$Date_Death)), n=nrow(data_outcome))
c(CFR_binomial$estimate,CFR_binomial$conf.int) |> c.text()

# Output unadjusted ('naive') CFR and 95% CI:
c.text(out_cfr$nCFR)

# Output delay adjusted CFR and 95% CI:
c.text(out_cfr$cCFR)


# Add comparison with EpiLine --------------------------------------------------

install_github("BDI-pathogens/EpiLine")
library(EpiLine)

# Filter on known onset and deaths
data_death <- data_ebola %>% filter(!is.na(Date_onset) & !is.na(Date_Death)) %>% arrange(Date_onset)

# Simple fitdistrplus fitting:
day_differences <- as.numeric(data_death$Date_Death-data_death$Date_onset)
param_index <- fitdistrplus::fitdist(day_differences,"lnorm") 

xx <- 1:40
hist(day_differences,prob=T,xlab="days",main="onset-to-death")
lines(xx,dlnorm(xx,param_index$estimate[["meanlog"]],param_index$estimate[["sdlog"]]))


# Define vectors
report_date <- min(data_outcome$Date_case,na.rm=T)
reported    <- case_data$count 
ll_symptom  <- as.numeric(data_death$Date_onset - report_date)
ll_outcome   <- as.numeric(data_death$Date_Death - report_date)

# fit using EpiLine model
mcmc_n_samples <- 100
mcmc_n_chains  <- 1
fit <- symptom_report.fit( reported, ll_symptom, ll_outcome, report_date = report_date, 
                           mcmc_n_samples = mcmc_n_samples, mcmc_n_chains = mcmc_n_chains )

# Plot posterior time series estimates
fit$plot.symptoms()

# Plot posterior parameters estimate
fit$plot.r() # Growth rate
fit$plot.symptom_report.dist() # Onset-to-death interval
