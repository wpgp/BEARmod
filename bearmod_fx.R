##### BEARmod v.0.81
#
# Basic Epidemic, Activity, and Response simulation model
#
# v0.81
# - Major speed increase in infection 
#
# v0.8
# - added capacity for multiple timesteps per day. TSinday defaults to 1, but defines the number of timesteps per day in the model (must be integer)
# - added capacity for a probability of moving per timestep (fitting to Vodafone data). Not accounted for if prob_move_per_TS = 0
#
# v0.7 updates:
# - Add percentage of exposed people who are infectious
#
# v0.65 updates:
# - Cleaned up inputs into model for easier pre-processing
# - Added capacity for time-dependent contact rates
#
# v0.6 updates:
# - fixed bug when recovery rate data are missing a patch for a day
# - patched bug that could lead to negative nInf values (however, the actual solution needs revisiting!)
#
#v0.5 updates:
# - Added functionality to input relative movement table
# - added functionality for time-variable recovery rates
#
# v.0.4 Updates: 
# - Model calibrated using HKU studies. 
# - Fixed transmission term to Poisson distribution 
# - Added "date" versatility to use non-contiguous dates
#
# This model runs a basic SEIR model, with stochastic exposure, incubation, recovery, and movement
# Disease spread occurs each day, based on the movement patterns from specific days from mobile phone-derived data.
#
#
#
# See run_model.R for working example.
#
# TO DO:
# - Add in infectious period for part of the exposed period (new category: exposed and infectious)
#
#
#
######

library(lubridate)

#This function creates the starting population
InitiatePop = function(pat_locator,initialInf,initialExp){
  NPat = dim(pat_locator)[1]
  list(
    nInitialInf = initialInf,
    nInitialExp = initialExp,
    nInf = initialInf,
    nExp = initialExp,
    nTotal = pat_locator$pop,
    names = pat_locator$patNames,
    IDs = pat_locator$patIDs,
    relativeInf = rep(1,NPat),
    nRecoveredToday = rep(0,NPat),
    nInfectedToday = rep(0,NPat),
    nExposedToday = rep(0,NPat),
    nInfMovedToday = rep(0,NPat),
    controlled = rep(0,NPat)
  )
}

##### Epidemic functions: exposure, infectivity, recovery  ####
recoveryTimeStep = function(HPop, recrate_values,current_day){
  recrate = subset(recrate_values,date == current_day)$recrate
  
 #print(recrate)#print(paste0("Day ",current_day, " recovery rate: ", recrate))
  for (i in 1:length(HPop$nInf)){
  HPop$nRecoveredToday[i]= sum(rbinom(HPop$nInf[i],1,recrate))
  HPop$nInf[i] = HPop$nInf[i] - HPop$nRecoveredToday[i]

  }
  #print(paste0("Number of people recovering: ",sum(HPop$nRecoveredToday)))
  HPop
}

exposedtoinfTimeStep = function(HPop, exp_to_infrate){
  #(exp_to_infrate)
  for (i in 1:length(HPop$nInf)){
    HPop$nInfectedToday[i]= sum(rbinom(HPop$nExp[i],1,exp_to_infrate))
    HPop$nInf[i] = HPop$nInf[i] + HPop$nInfectedToday[i]
    HPop$nExp[i] = HPop$nExp[i] - HPop$nInfectedToday[i]
  }
  #print(paste0("Number of people newly infectious: ",sum(HPop$nInfectedToday)))
  HPop
}

exposedTimeStep = function(HPop, exposerate_df, current_day, exposed_pop_inf_prop){
  
  if (is.numeric(exposerate_df)){
    exposerate = exposerate_df
  }
  if (is.data.frame(exposerate_df)){
    exposerate = subset(exposerate_df, date == current_day)$exposerate
  }
  for (i in 1:length(HPop$nInf)){
    infectious_pop = HPop$nInf[i] + exposed_pop_inf_prop * HPop$nExp[i]
    #HPop$nExposedToday[i]= sum(rbinom(infectious_pop,1,exposerate)) * (1 - ( (HPop$nInf[i] + HPop$nExp[i]) / HPop$nTotal[i]))
    HPop$nExposedToday[i]= sum(rpois(infectious_pop,exposerate)) * (1 - ( (HPop$nInf[i] + HPop$nExp[i]) / HPop$nTotal[i]))


    if (HPop$nExp[i] + HPop$nExposedToday[i] < HPop$nTotal[i] - HPop$nInf[i]) {
      HPop$nExp[i] = HPop$nExp[i] + HPop$nExposedToday[i]
      
    } else {
      HPop$nExposedToday[i] = HPop$nTotal[i] - HPop$nInf[i] - HPop$nExp[i] 
      HPop$nExp[i]  = HPop$nTotal[i] - HPop$nInf[i]
      
    }
  }
  #print(paste0("Number of people newly exposed: ",sum(HPop$nExposedToday)))
  HPop
}




####### Activity functions: Human movement ####
movementTimeStep = function(HPop, mobmat,day,control_df,prob_move_per_TS){
  movement_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  daily_move = subset(mobmat,date == day)
  
  daily_move = subset(daily_move,!is.na(fr_pat) & !is.na(to_pat) & !is.na(fr_users) & !is.na(movers))
  
  daily_move_mat = daily_move[,is.element(names(daily_move),c("fr_pat","to_pat","fr_users","movers"))]
  daily_move_mat = as.matrix(daily_move_mat)
  col1 = which(colnames(daily_move_mat) == "fr_pat")
  col2=which(colnames(daily_move_mat) == "to_pat")
  
  colmove = which(colnames(daily_move_mat) == "movers")
  colusers=which(colnames(daily_move_mat) == "fr_users")
  movement_matrix[daily_move_mat[,c(col1,col2)]] = daily_move_mat[,colmove]/daily_move_mat[,colusers]
  if (length(which(rowSums(movement_matrix)>1)) > 0){
    print("Warning: row sums > 1 in movement matrix. Correcting...")
    correctingrows = which(rowSums(movement_matrix)>1)
    for (i in correctingrows){
    movement_matrix[i,] = movement_matrix[i,] /sum(movement_matrix[i,] )
    }
  }
  if (prob_move_per_TS > 0){
    movement_matrix = movement_matrix*prob_move_per_TS
  }
  for (i in 1:length(patIDs)){
    movement_matrix[i,i] = 1 - sum(movement_matrix[i,-i])
  }
  
  HPop$controlled = rep(0,length(HPop$names))
  if (length(which(control_df$date == day)) > 0){
    control_df_sub = subset(control_df,date == day)
    if (dim(control_df_sub)[1] > 0){
    for (i in 1:dim(control_df_sub)[1]){
      HPop$controlled[which(HPop$names == control_df_sub$from[i])] = control_df_sub$relative_move[i]
      
    }
  }
  }
  if (sum(HPop$controlled)>0){
    movement_matrix = stopMovement(HPop,movement_matrix,day)
  }
  #deterministic version
  #HPop$nInfMovedToday = colSums(diag(HPop$nInf) %*% movement_matrix) - HPop$nInf
  #HPop$nInf = colSums(diag(HPop$nInf) %*% movement_matrix)
  HPop$nInf = ceiling(HPop$nInf)
  # stochastic version
    z <- rbinom(n=NPat^2,size = rep(HPop$nInf,each=NPat),prob = t(movement_matrix)[])
  moved_matrix = t(matrix(z,NPat,NPat,dimnames=list(patIDs,patIDs)))
  for (i in 1:dim(moved_matrix)[1]){
     if (sum(moved_matrix[i,]) > 0){
     moved_matrix[i,] = moved_matrix[i,]/sum(moved_matrix[i,]) * HPop$nInf[i]
     }
  }
  HPop$nInfMovedToday = ceiling(colSums(moved_matrix))
  
  HPop$nInf = HPop$nInf - floor(rowSums(moved_matrix)) + ceiling(colSums(moved_matrix))
  #quick fix
  for (i in 1:length(HPop$nInf)){
    if (HPop$nInf[i] > HPop$nTotal[i]){
      HPop$nInf[i] = HPop$nTotal[i]
    }
  }
  
  #print(paste0("Number of infected people moving: ",sum(abs(HPop$nInfMovedToday))/2))
  HPop
}


###### Response functions: Control
#relative_movement is the proportion of original movement out/in that we want to keep -- ie. .1 = 10% of original movement rate
stopMovement = function(HPop,mobmat,current_date){
  stopping = which(HPop$controlled > 0)
    if (length(stopping) > 0){
     # print(paste("stopping movement in patches", HPop$names[stopping]))
      for (ctrl_pat in stopping){
    control_patches = HPop$IDs[ctrl_pat]
    mobmat[control_patches,] = mobmat[control_patches,] * HPop$controlled[ctrl_pat]
    mobmat[,control_patches] = mobmat[,control_patches] * HPop$controlled[ctrl_pat]
    for (i in 1:length(HPop$IDs)){
      mobmat[i,i] = 1 - sum(mobmat[i,-i])
    }
      }
    }
    mobmat
}




###### Master function  ####
runSim = function(HPop,pat_info,control_info,mobmat,day_list,recrate_values,exposerate_df,exposepd,exposed_pop_inf_prop = 0,TSinday = 1,prob_move_per_TS=0) {
  
  
  epidemic_curve <- data.frame(Date=as.Date(character()),
                               inf=c(),
                               stringsAsFactors=FALSE) 
  
  if (TSinday > 1){
    #recrate_values$recrate = 1-(1-recrate_values$recrate)^(1/TSinday)
    exposetoinfrate = 1/exposepd
    exposepd = 1/(1 - exp(log(1-exposetoinfrate) / TSinday))
    #recrate_values$recrate = 1 - ((1 - recrate_values$recrate) ^ (1/TSinday))
    recrate_values$recrate = 1 - exp(log(1-recrate_values$recrate) / TSinday)
    if (is.numeric(exposerate_df)){
     # exposerate_df = 1-(1-exposerate_df)^(1/TSinday)
      exposerate_df = exposerate_df/TSinday      
     # recrate_values$recrate = 1 - ((1 - recrate_values$recrate) ^ (1/TSinday))
    }
    if (is.data.frame(exposerate_df)){
     # exposerate_df$exposerate = 1-(1-exposerate_df$exposerate)^(1/TSinday)
      exposerate_df$exposerate = 1 - exp(log(1-exposerate_df$exposerate) / TSinday)
    }
  }

  all_spread = matrix(0,length(day_list),length(HPop$nInf))
  colnames(all_spread) = HPop$names
  #print(all_dates)
  for (current_day in 1:length(day_list)){
    for (current_TS in 1:TSinday){
    print(day_list[current_day])
    HPop = recoveryTimeStep(HPop,recrate_values,day_list[current_day])
    HPop = exposedtoinfTimeStep(HPop,1/exposepd)
    HPop = exposedTimeStep(HPop,exposerate_df, day_list[current_day], exposed_pop_inf_prop)
    
    HPop = movementTimeStep(HPop,mobmat,day_list[current_day],control_info,prob_move_per_TS)
    }
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = day_list[current_day], inf = sum(HPop$nInf)))
    all_spread[current_day,] = HPop$nInf
    
  }
  all_spread_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_2= cbind(all_spread_2,all_spread)
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread_2)
}


