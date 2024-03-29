# - - - - - - - - - - - - - - - - - - - - - - - 
# Deconvolution of simulated infection data
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load libraries
library(MASS)
library(tidyverse)
install_github("epiverse-trace/epiparameter")
library(epiparameter)

# Simulate infection dynamics
xx <- 0:150
data_infections <- 0.1*(2+sin(8*pi*(xx-20)/365))
n_inf <- length(data_infections) # number of days to consider
data_infections <- data_infections #* rlnorm(n_inf,0,0.2) # add some noise

# Set delay function pmf
#p_by_day <- #epiparameter::epidist("SARS_CoV_2_wildtype","incubation")$pmf
mean_p <- 5
scale_p <- 1
shift_p <- 0

p_by_day <- function(x){dgamma(x,shape=mean_p/scale_p,scale=scale_p)}
plot(1:20,p_by_day(1:20))

# Define transition matrix to construct outcome data
f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
n_delay_days <- 50 # maximum delay period to consider

for(ii in 1:n_inf){
  i_max <- min(ii+n_delay_days-1,n_inf)
  j_max <- min(n_inf-ii+1,n_delay_days)
  
  f_matrix[ii:i_max,ii] <- p_by_day(0:(j_max-1)) # fill matrix entries
  
}


# Quick simulation --------------------------------------------------------

# Simulate outcomes
data_outcomes <- f_matrix %*% data_infections

par(mfrow=c(1,1),mgp=c(2,0.7,0),mar = c(3,3,1,1))

# Plot original incidence
plot(data_infections,yaxs="i",ylab="daily incidence (%)",ylim=c(0,0.5),xlab="days")
points(data_outcomes,col="red")
title(main=LETTERS[1],adj=0);letter_x <- letter_x+1



# Run deconvolution methods -----------------------------------------------


# Inversion function
invert_f <- ginv(f_matrix)

# Function to estimate infection incidence from delayed outcomes
estimate_infections <- function(delayed_outcomes){
  
  # Define transition matrix - 
  n_inf <- length(delayed_outcomes)
  f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
  n_delay_days <- 30 # maximum incubation period to consider
  
  for(ii in 1:n_inf){
    i_max <- min(ii+n_delay_days-1,n_inf)
    j_max <- min(n_inf-ii+1,n_delay_days)
    
    f_matrix[ii:i_max,ii] <- p_by_day(0:(j_max-1))
    
  }
  
  # Calculate Moore-Penrose pseudoinverse matrix 
  invert_f <- ginv(f_matrix)
  
  # Apply inversion matrix
  output <- invert_f %*% delayed_outcomes
  output
  
}

# Run simulation function to generate delayed outcomes from original infections
data_outcomes <- f_matrix %*% data_infections
#data_outcomes <- data_outcomes * rlnorm(length(data_outcomes),0,0.05) # add some noise

par(mfrow=c(1,2),mgp=c(2,0.7,0),mar = c(3,3,1,1))
letter_x <- 1

# Plot original incidence
plot(data_infections,yaxs="i",ylab="daily incidence (%)",ylim=c(0,0.5),xlab="days")
points(data_outcomes,col="red")
title(main=LETTERS[1],adj=0);letter_x <- letter_x+1

inc1 <- estimate_infections(data_outcomes)
lines(inc1,col="blue",lwd=2)

x_shift <- 60; y_shift <- 0.45
text(labels="true infections",x=x_shift,y=y_shift,adj=0,col="black")
text(labels="delayed outcomes",x=x_shift,y=(y_shift-0.025),adj=0,col="red")
text(labels="estimated infections",x=x_shift,y=(y_shift-0.05),adj=0,col="blue")

# Run again but with noise on delayed outcome observations

data_outcomes2 <- data_outcomes * rlnorm(length(data_outcomes),0,0.01) # add some noise


# Plot original incidence
plot(data_infections,yaxs="i",ylab="daily incidence (%)",ylim=c(0,0.5),xlab="days")
points(data_outcomes2,col="red")
title(main=LETTERS[2],adj=0);letter_x <- letter_x+1

inc2 <- estimate_infections(data_outcomes2)
lines(inc2,col="blue",lwd=2)
