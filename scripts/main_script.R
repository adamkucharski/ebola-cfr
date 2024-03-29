# - - - - - - - - - - - - - - - - - - - - - - - 
# CFR estimation analysis
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load libraries
library(tidyverse)

# Define probability mass function for onset-to-death
# Source: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)31387-4/fulltext
onset_to_death_ebola <- function(x){dgamma(x,shape=1.651,scale=1/0.202)}

# Load data and functions
source("R/data_load.R")
source("R/cfr_function.R")

# Run example delay adjusted CFR estimation on 1976 Ebola data: -----------

data_in <- data1976 # Incidence data
correct_cfr <- correct_cfr_1976 # Actual final CFR

# Estimate CFR over first 30 days
cutoff <- 30
out_cfr <- scale_cfr(case_incidence = data_in$cases[1:cutoff],
                     death_incidence = data_in$deaths[1:cutoff],
                     delay_fun = onset_to_death_ebola)

# Unadjusted ('naive') CFR and 95% CI:
out_cfr$nCFR

# Delay adjusted CFR and 95% CI:
out_cfr$cCFR


# Plot estimates over time on both 1976 and 2014 data ------------------------------------------

# Black lines, daily new cases; orange lines, daily new deaths
# Red lines, naive CFR calculation, i.e. deaths/cases, with 95% CI; blue lines, delay-adjusted calculation with 95%

source("R/plot_examples.R")