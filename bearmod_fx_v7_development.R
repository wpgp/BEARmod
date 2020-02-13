##### BEARmod v.0.6
#
# Basic Epidemic, Activity, and Response simulation model
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
# Disease spread occurs each day, based on the movement patterns from specific days from Baidu data.
#
#
#
# See run_model_example.R for fully working example.
#
# TO DO:
# - Add in infectious period for part of the exposed period (new category: exposed and infectious)
# - Add in ability to reduce contact rates in certain places
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
    nTotal = pat_locator$TOTAL_POP2015,
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
  #print(paste0("Day ",current_day, " recovery rate: ", recrate))
  for (i in 1:length(HPop$nInf)){
  HPop$nRecoveredToday[i]= sum(rbinom(HPop$nInf[i],1,recrate))
  HPop$nInf[i] = HPop$nInf[i] - HPop$nRecoveredToday[i]

  }
  #print(paste0("Number of people recovering: ",sum(HPop$nRecoveredToday)))
  HPop
}

exposedtoinfTimeStep = function(HPop, exp_to_infrate){
  
  for (i in 1:length(HPop$nInf)){
    HPop$nInfectedToday[i]= sum(rbinom(HPop$nExp[i],1,exp_to_infrate))
    HPop$nInf[i] = HPop$nInf[i] + HPop$nInfectedToday[i]
    HPop$nExp[i] = HPop$nExp[i] - HPop$nInfectedToday[i]
  }
  #print(paste0("Number of people newly infectious: ",sum(HPop$nInfectedToday)))
  HPop
}

exposedTimeStep = function(HPop, exposerate){
  for (i in 1:length(HPop$nInf)){
    HPop$nExposedToday[i]= sum(rpois(HPop$nInf[i],exposerate)) * (1 - ( (HPop$nInf[i] + HPop$nExp[i]) / HPop$nTotal[i]))
  if ( HPop$nExp[i] + HPop$nExposedToday[i] < HPop$nTotal[i] - HPop$nInf[i]) {
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
movementTimeStep = function(HPop, mobmat,day,control_df){
  movement_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  daily_move = subset(mobmat,Date == day)
  daily_move_mat = as.matrix(daily_move[,c(4:8)])
  movement_matrix[daily_move_mat[,c(2,3)]] = daily_move_mat[,1]/daily_move_mat[,4]
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
  moved_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  for (i in 1:dim(moved_matrix)[1]){
    #for (j in 1:dim(moved_matrix)[2]){
    moved_matrix[i,] = rbinom(length(moved_matrix[i,]),HPop$nInf[i],prob = movement_matrix[i,])
    if (sum(moved_matrix[i,]) > 0){
    moved_matrix[i,] = moved_matrix[i,]/sum(moved_matrix[i,]) * HPop$nInf[i]
    }

  }
  
  HPop$nInfMovedToday = ceiling(colSums(moved_matrix))
  HPop$nInf = HPop$nInf - floor(rowSums(moved_matrix)) + ceiling(colSums(moved_matrix))
  
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
runSim = function(HPop,pat_info,control_info,mobmat,day_list,recrate_values,exposerate,exposepd) {
  
  
  epidemic_curve <- data.frame(Date=as.Date(character()),
                               inf=c(),
                               stringsAsFactors=FALSE) 
  

  all_spread = matrix(0,length(day_list),length(HPop$nInf))
  colnames(all_spread) = HPop$names
  #print(all_dates)
  for (current_day in 1:length(day_list)){
    
   # print(day_list[current_day])
    HPop = recoveryTimeStep(HPop,recrate_values,day_list[current_day])
    HPop = exposedtoinfTimeStep(HPop,1/exposepd)
    HPop = exposedTimeStep(HPop,exposerate)
    HPop = movementTimeStep(HPop,mobmat,day_list[current_day],control_info)
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = day_list[current_day], inf = sum(HPop$nInf)))
    all_spread[current_day,] = HPop$nInf
    
  }
  all_spread_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_2= cbind(all_spread_2,all_spread)
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread_2)
}


