# - - - - - - - - - - - - - - - - - - - - - - - 
# Load previous Ebola outbreak data
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load data
data1976 <- read_csv("data/ebola_1976.csv")
correct_cfr_1976 <- 234/245 # calculate final CFR for validation

data2014 <- read_csv("data/ebola_2014.csv")
correct_cfr_2014 <- 0.7

# Convert cumulative data into estimates of incidence
cases_d_2014 <- rep(diff(data2014$cases)/as.numeric(diff(data2014$date)),as.numeric(diff(data2014$date))) |> round()
deaths_d_2014 <- rep(diff(data2014$deaths)/as.numeric(diff(data2014$date)),as.numeric(diff(data2014$date))) |> round()
dates_d_2014 <- min(data2014$date) + 0:(length(cases_d_2014)-1)

# Compile daily data
data2014_daily <- tibble(date=dates_d_2014,cases=cases_d_2014,deaths=deaths_d_2014)

# Compile data in coarseDataTools format:
data1976_all <- read_csv("data/ebola_1976_all_outcomes.csv")
data1976_all <- data1976_all |> mutate(new.times = as.numeric(time - min(data1976_all$time)))
