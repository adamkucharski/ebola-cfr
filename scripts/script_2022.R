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

# Define probability mass function for onset-to-death
onset_to_death_ebola <- epiparameter::epidist("ebola","onset_to_death")$pmf

xx <- 0:20
plot(xx,onset_to_death_ebola(xx),xlab="days")

# Load data and functions
source("R/cfr_function.R")


# Import 2022 Ebola data --------------------------------------------------
data_ebola <- read_csv("https://raw.githubusercontent.com/globaldothealth/ebola/main/latest.csv")

# Estimate CFR for cases with known outcomes --------------------------------------------------
data_ebola <- data_ebola |> mutate(Date_case = as.Date(ifelse(!is.na(Date_onset),Date_onset,Date_confirmation),"1970-01-01"))
data_outcome <- data_ebola %>% filter(!is.na(Date_death) | !is.na(Date_recovered))

# Output estimate based on known outcomes
CFR_binomial <- prop.test(x=sum(!is.na(data_outcome$Date_death)), n=nrow(data_outcome))
c(CFR_binomial$estimate,CFR_binomial$conf.int) |> c.text()

# Estimate CFR from incidence data --------------------------------------------------
case_data <- incidence2::incidence(data_ebola,Date_case) |> complete_counts()
death_date <- incidence2::incidence(data_ebola,Date_death) |> complete_counts()

# Estimate CFR:
out_cfr <- scale_cfr(case_incidence = case_data$count,
                     death_incidence = death_date$count,
                     delay_fun = onset_to_death_ebola)

# Output unadjusted ('naive') CFR and 95% CI:
c.text(out_cfr$nCFR)

# Output delay adjusted CFR and 95% CI:
c.text(out_cfr$cCFR)


