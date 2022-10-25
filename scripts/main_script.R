# Load libraries and data
library(tidyverse)

# Define probability mass function for onset-to-death
# Source: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)31387-4/fulltext
onset_to_death_ebola <- function(x){dgamma(x,shape=1.651,scale=1/0.202)}

# Load functions
source("R/data_load.R")
source("R/cfr_function.R")
source("R/plot_examples.R")

# Run example correction on 1976 Ebola data:
data_in <- data1976 # Incidence data
correct_cfr <- correct_cfr_1976 # Actual final CFR

dates_d <- data_in$date
cases_d <- data_in$cases
deaths_d <- data_in$deaths

# Run on all data

out_cfr <- scale_cfr(case_incidence = cases_d,
                     death_incidence = deaths_d,
                     delay_fun = onset_to_death_ebola)

out_cfr