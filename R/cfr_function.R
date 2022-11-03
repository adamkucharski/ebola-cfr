# - - - - - - - - - - - - - - - - - - - - - - - 
# CFR MLE estimation function
# Author: Adam Kucharski
# https://github.com/adamkucharski/ebola-cfr
# - - - - - - - - - - - - - - - - - - - - - - - 

scale_cfr <- function(case_incidence, # daily case incidence
                      death_incidence, #daily death incidence
                      delay_fun # onset-to-death distribution pmf
                      ){
  
  cumulative_known_t <- 0 # store cumulative cases with known outcome at time tt
  
  # Sum over cases up to time tt
  for(ii in 1:length(case_incidence)){
    known_i <- 0 # store number of cases with known clinical outcome at time ii
    for(jj in 0:(ii - 1)){ # iterate over cases up to this point
      known_jj <- (case_incidence[ii - jj]*delay_fun(jj))
      known_i <- known_i + known_jj
    }
    cumulative_known_t <- cumulative_known_t + known_i # tally cumulative known outcomes
  }

  D_t <- sum(death_incidence) # cumulative deaths
  C_t <- sum(case_incidence) # cumulative cases
  
  u_t <- cumulative_known_t/sum(case_incidence) # proportion of cases with known outcome
  
  # - - - 
  # Calculate naive CFR value and 95% binomial interval
  htest <- binom.test(D_t,C_t, p = 1,conf.level=0.95)
  b_t <- c(D_t/C_t,htest$conf.int[1:2])
  
  # - - - 
  # MLE estimation for corrected CFR
  pprange <- seq(1e-3,1,1e-3)
  lik <- matrix(NA,nrow=length(pprange))
  
  for(ii in 1:length(pprange)){
    p_t <- pprange[ii]
    
    # Calculate likelihood - use binomial for small samples and Poisson approximation for larger numbers
    if(C_t<200){
      lik[ii] <- log(choose(round(u_t*C_t),D_t))+ D_t*log(p_t)+(u_t*C_t-D_t)*log(1-p_t)
    }else{
      lik[ii] <- dpois(D_t,p_t*u_t*C_t,log=T)
    }
    
  }
  
  mid_val <- pprange[lik==max(lik)] # MLE estimate
  ci_95 <- pprange[lik>=(max(lik)-1.92)] 
  ci_95 <- c(min(ci_95),max(ci_95)) # 95% interval
  
  if(is.na(max(lik))){mid_val <- NA;ci_95 <- c(NA,NA)}
  
  data.frame(nCFR = b_t, cCFR = c(mid_val,ci_95), total_deaths = D_t, 
             total_cases = C_t,proportion_cases_known = u_t)
  
  
}


# Output estimate and 95% as text -----------------------------------------

c.text <- function(x,sigF=3){
  y <- signif(x,sigF)
  paste(y[1]," (",y[2],"-",y[3],")",sep="")
}
