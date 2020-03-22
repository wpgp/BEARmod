####Model running code for BEARmod v.0.6
rm(list=ls())
library(data.table) # fread - fastly reading data
library(lubridate)

source("bearmod_fx.R")
  source("preprocess_small.R")
#Initial parameters
NPat = length(patNames)
patnInf = rep(0,NPat)
patnExp = c(rep(0,NPat) )

pat_locator$pop = 100

#start infection in Wuhan
patnInf[which(patNames == 1)] = 50
#recovery rate variable
recover_df = data.frame(date = seq(from=min(movement_data$date),to=max(movement_data$date),by="days"),recrate = recrate)
relative_move_data=data.frame()
 
#### Running the model  ####



HPop = InitiatePop(pat_locator,patnInf,patnExp)
###dates of simulation


 input_dates = rep("2020-05-01",50)
results = list()

for (run in 1:500){
  
  HPop_update = runSim(HPop,pat_locator,relative_move_data,movement_data, input_dates,recover_df, exposerate,exposepd,exposed_pop_inf_prop = 0, TSinday = 1)
  print(paste0("Run # ",run))
  results[[run]] = HPop_update$all_spread
}
save(results,file="results.RData")
# 