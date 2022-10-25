# Plot results of 1976 and 2014 Ebola epidemics ---------------------------

par(mfcol=c(2,2),mar=c(4,4,1,1),las=1)

for(ii in 1:2){
  
  if(ii==1){data_in <- data1976; correct_cfr <- correct_cfr_1976; labelx <- "Ebola 1976"}
  if(ii==2){data_in <- data2014_daily; correct_cfr <- correct_cfr_2014; labelx <- "Ebola 2014"}
  
  dates_d <- data_in$date
  cases_d <- data_in$cases
  deaths_d <- data_in$deaths
  
  # Run on all data

  # Show rolling window for data
  merge_data <- data.frame(cases_d,deaths_d)
  burn_in <- round(0.15*length(cases_d)) # Omit initial data points
  time_data <- burn_in:nrow(merge_data) # Focus from day=burn_in onwards
  
  #time_data <- 11:35
  out_all <- sapply(time_data,function(x){
    out_cfr <- scale_cfr(case_incidence = merge_data[1:x,]$cases_d,
                         death_incidence = merge_data[1:x,]$deaths_d,
                         delay_fun = onset_to_death_ebola)
    c(out_cfr$nCFR,out_cfr$cCFR)
  })
  
  out_all <- data.frame(t(out_all)); names(out_all) <- c("mid_n","lower_n","upper_n","mid_c","lower_c","upper_c")
  
  # Plot comparative results for
  dates_burn_in <- tail(dates_d,-burn_in+1)
  date_lim <- c(min(dates_d),max(dates_d))
  
  plot(dates_d,cases_d,xlim=date_lim,ylim=c(0,1.2*max(cases_d)),
       type="s",lwd=2,col="black",ylab="cases and deaths",xlab="",yaxs="i",main=labelx) # Plot naive CFR
  lines(dates_d,deaths_d,type="s",lwd=2,col="orange") # Plot naive CFR
  
  plot(dates_burn_in,out_all$mid_n,xlim=date_lim,ylim=c(0,1),xlab="",yaxs="i",type="l",lwd=2,col="red",ylab="CFR") # Plot naive CFR
  lines(dates_burn_in,out_all$lower_n,col="red");lines(dates_burn_in,out_all$upper_n,col="red",lwd=0.5)
  
  lines(dates_burn_in,out_all$mid_c,col="blue",lwd=2) # Plot corrected CFR
  lines(dates_burn_in,out_all$lower_c,col="blue");lines(dates_burn_in,out_all$upper_c,col="blue",lwd=0.5)
  
  lines(min(dates_d)+c(0,1e3),c(correct_cfr,correct_cfr),lty=2) # Plot final CFR
  
}

dev.copy(pdf,"plots/CFR_figure.pdf",width=10,height=7) 
dev.off()