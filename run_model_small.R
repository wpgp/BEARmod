####Model running code for BEARmod v.0.6
rm(list=ls())
library(data.table) # fread - fastly reading data
library(lubridate)

# setwd('//worldpop.files.soton.ac.uk/Worldpop/Projects/WP519091_Seasonality')
# setwd('D:/OneDrive - University of Southampton/Wuhan Coronavirus R0/Spread risk')
#setwd('C:/Users/sl4m18/OneDrive - University of Southampton/Wuhan Coronavirus R0/Spread risk')


source("bearmod/BEARmod_development/bearmod_fx_dev.R")
# source("bearmod/bearmod_fx.R")
source("bearmod/BEARmod_development/preprocess_small.R")
#Initial parameters
NPat = length(patNames)
patnInf = rep(0,NPat)
patnExp = c(rep(0,NPat) )

pat_locator$pop = 100

#start infection in Wuhan
patnInf[which(patNames == 1)] = 50
#recovery rate variable
recover_df = data.frame(date = seq(from=min(movement_data$date),to=max(movement_data$date),by="days"),recrate = recrate)
relative_move_df=data.frame()
 
#### Running the model  ####



HPop = InitiatePop(pat_locator,patnInf,patnExp)
###dates of simulation


 input_dates = rep("2020-05-01",50)
 # input_dates = seq(date("2013-12-02"),date("2014-2-13"),by="days") # coresponding to the period from 2020-12-08 to 2 wks after LNY's day 
# input_dates = seq(date("2013-12-02"),date("2014-2-27"),by="days") # coresponding to the period from 2020-12-08 to 4 wks after LNY's day
results = list()

for (run in 1:500){
  
  HPop_update2 = runSim(HPop,pat_locator,relative_move_data,movement_data, input_dates,recover_df, exposerate,exposepd,exposed_pop_inf_prop = 0, TSinday = 1)
  print(paste0("Run # ",run))
  results[[run]] = HPop_update$all_spread
}
save(results,file="results.RData")
# 