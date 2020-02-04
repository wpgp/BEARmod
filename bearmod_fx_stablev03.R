##### BEARmod v.0.3
#
# Basic Epidemic, Activity, and Response simulation model
#
# This model runs a basic SEIR model, with stochastic exposure, incubation, recovery, and movement
# Disease spread occurs each day, based on the movement patterns from specific days from Baidu data
#
# TO DO:
# - Add in infectious period for part of the exposed period (new category: exposed and infectious)
# - Add in ability to reduce contact rates in certain places
#
#
#
# - Should results (# people infected, etc) be integers?
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
    nInfMovedToday = rep(0,NPat)
  )
}

##### Epidemic functions: exposure, infectivity, recovery  ####
recoveryTimeStep = function(HPop, recrate){
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
    HPop$nExposedToday[i]= round(sum(rbinom(HPop$nInf[i],1,exposerate)) * (1 - ( (HPop$nInf[i] + HPop$nExp[i]) / HPop$nTotal[i])))
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
movementTimeStep = function(HPop, mobmat,day,control_df,relative_movement){
  movement_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  daily_move = subset(mobmat,Date == day)
  daily_move_mat = as.matrix(daily_move[,c(4:8)])
  movement_matrix[daily_move_mat[,c(2,3)]] = daily_move_mat[,1]/daily_move_mat[,4]
  for (i in 1:length(patIDs)){
    movement_matrix[i,i] = 1 - sum(movement_matrix[i,-i])
  }
  if (length(which(control_df$control_date < day)) > 0){
    movement_matrix = stopMovement(control_df,movement_matrix,day,relative_movement)
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
  
  HPop$nInfMovedToday = colSums(moved_matrix)
  HPop$nInf = HPop$nInf - rowSums(moved_matrix) + colSums(moved_matrix)
  
  #print(paste0("Number of infected people moving: ",sum(abs(HPop$nInfMovedToday))/2))
  HPop
}


###### Response functions: Control
#relative_movement is the proportion of original movement out/in that we want to keep -- ie. .1 = 10% of original movement rate
stopMovement = function(control_df,mobmat,current_date,relative_movement){
  stopping = which(control_df$control_date < current_date)
    if (length(stopping) > 0){
    #  print(paste("stopping movement in patches", control_df$patIDs[stopping]))
    control_patches = control_df$patIDs[stopping]
    mobmat[control_patches,] = mobmat[control_patches,] * relative_movement
    mobmat[,control_patches] = mobmat[,control_patches] * relative_movement
    for (i in 1:length(patIDs)){
      mobmat[i,i] = 1 - sum(mobmat[i,-i])
    }
    }
    mobmat
}




###### Master function  ####
runSim = function(HPop,pat_info,mobmat,day_start,day_end,recrate,exposerate,exposepd,relative_movement = 1) {
  
  
  epidemic_curve <- data.frame(Date=as.Date(character()),
                               inf=c(),
                               stringsAsFactors=FALSE) 
  
  all_dates = seq(day_start,day_end,by="days")
  
  all_spread = matrix(0,length(HPop$nInf),length(all_dates))
  #print(all_dates)
  for (current_day in 1:length(all_dates)){
    
    #print(all_dates[current_day])
    HPop = recoveryTimeStep(HPop,recrate)
    HPop = exposedtoinfTimeStep(HPop,1/exposepd)
    HPop = exposedTimeStep(HPop,exposerate)
    HPop = movementTimeStep(HPop,mobmat,all_dates[current_day],pat_info,relative_movement)
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = all_dates[current_day], inf = sum(HPop$nInf)))
    all_spread[,current_day] = HPop$nInf
    
  }
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread)
}


