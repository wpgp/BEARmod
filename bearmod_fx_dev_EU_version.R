##### BEARmod v.0.81
#
# Basic Epidemic, Activity, and Response simulation model
#
#
# v0.92
# #fixed rounding error bugs throughout
#
#
# v0.91
# - fixed issue with movement that could lead to more people moving into a patch than there are total. for now, this was fixed by making "room" for infected people by removing recovered people.
# 
#
# v0.9
# - added functionality for time spent version of movement matrix (for now implemented in a new function, runSim_timespent())
# - fixed bug in how model tracks recovered people
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
    nRec = rep(0,NPat),
    nTotal = pat_locator$pop,
    names = as.character(pat_locator$patNames),
    IDs = pat_locator$patIDs,
    relativeInf = rep(1,NPat),
    nRecoveredToday = rep(0,NPat),
    nInfectedToday = rep(0,NPat),
    nExposedToday = rep(0,NPat),
    nInfMovedToday = rep(0,NPat),
    relativeContactRate = rep(1,NPat),
    nImportedTodayIntl = rep(0,NPat),
    nExportedTodayIntl = rep(0,NPat),
    nImportedToday = rep(0,NPat),
    nExportedToday = rep(0,NPat),
    controlled = rep(0,NPat)
  )
}

##### Epidemic functions: exposure, infectivity, recovery  ####
recoveryTimeStep = function(HPop, recrate_values,current_day){
  recrate = subset(recrate_values,date == current_day)$recrate

  HPop$nInf = round(HPop$nInf)
 #print(recrate)#print(paste0("Day ",current_day, " recovery rate: ", recrate))
  for (i in 1:length(HPop$nInf)){
    #print(HPop$nInf[i])
   # print(recrate)
    recrate_new = recrate
    #if a nuts3 area had fewer than 100 active cases, we assumed test-and-trace would be able to double the recovery rate (helped quicken pace to local elimination in our simulations)
   if (HPop$nInf[i] < 100){
     if (current_day > "2020-03-30"){
        recrate_new = recrate_new * 2
     }
    }
  HPop$nRecoveredToday[i]= rbinom(1,HPop$nInf[i],recrate_new)
#print(HPop$nRecoveredToday[i])
  HPop$nInf[i] = HPop$nInf[i] - HPop$nRecoveredToday[i]
  HPop$nRec[i] = HPop$nRec[i] + HPop$nRecoveredToday[i]
  
  }
  for (i in 1:length(HPop$nInf)){
    if (HPop$nRec[i] > HPop$nTotal[i] - HPop$nExp[i] - HPop$nInf[i]){
      HPop$nInf[i] = max(0,HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec[i])
      HPop$nRec[i] = max(0, HPop$nTotal[i] - HPop$nExp[i]- HPop$nInf[i])
      
    }
    if (HPop$nRec[i] <0 ){
      
      HPop$nRec[i] = 0
    }
  }

  #print(paste0("Number of people recovering: ",sum(HPop$nRecoveredToday)))
  HPop
}

exposedtoinfTimeStep = function(HPop, exp_to_infrate){
  #(exp_to_infrate)
  HPop$nExp = round(HPop$nExp)
  for (i in 1:length(HPop$nInf)){
    #print(HPop$nExposedToday[i])
    HPop$nInfectedToday[i]= sum(rbinom(HPop$nExp[i],1,exp_to_infrate))
    #if (HPop$nInf[i] + HPop$nInfectedToday[i] < HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec[i] ) {
      HPop$nInf[i] = HPop$nInf[i] + HPop$nInfectedToday[i]
      
   # } else {
    #  HPop$nInfectedToday[i] = max(0,HPop$nTotal[i] - HPop$nInf[i] - HPop$nExp[i]- HPop$nRec[i])
      
    #  HPop$nInf[i]  = max(0,HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec[i])
      
   # }
    HPop$nExp[i] = HPop$nExp[i] - HPop$nInfectedToday[i]
    
  }
  #print(paste0("Number of people newly infectious: ",sum(HPop$nInfectedToday)))
  HPop
}

exposedTimeStep = function(HPop, exposerate_df,current_day, exposed_pop_inf_prop){
  
  if (is.numeric(exposerate_df)){
    exposerate = exposerate_df
  }
  if (is.data.frame(exposerate_df)){
    exposerate = subset(exposerate_df, date == current_day)$exposerate
  }
  exposerate_relative = HPop$relativeContactRate * exposerate
  print(mean(HPop$relativeContactRate))
  for (i in 1:length(HPop$nInf)){
    infectious_pop = HPop$nInf[i] + exposed_pop_inf_prop * HPop$nExp[i]
    #print(infectious_pop)
    #HPop$nExposedToday[i]= sum(rbinom(infectious_pop,1,exposerate)) * (1 - ( (HPop$nInf[i] + HPop$nExp[i]) / HPop$nTotal[i]))
    HPop$nExposedToday[i]= sum(rpois(infectious_pop,exposerate_relative[i])) * (1 - min(1, ( (HPop$nInf[i] + HPop$nExp[i] + HPop$nRec[i]) / HPop$nTotal[i]) ))
    
    
    if (HPop$nExp[i] + HPop$nExposedToday[i] < HPop$nTotal[i] - HPop$nInf[i] - HPop$nRec[i] ) {
      HPop$nExp[i] = HPop$nExp[i] + HPop$nExposedToday[i]
    
    } else {
      HPop$nExposedToday[i] = max(0,HPop$nTotal[i] - HPop$nInf[i] - HPop$nExp[i]- HPop$nRec[i])
      
      HPop$nExp[i]  = max(0,HPop$nTotal[i] - HPop$nInf[i] - HPop$nRec[i])
      
    }
  }
  #print(paste0("Number of people newly exposed: ",sum(HPop$nExposedToday)))
  HPop
}





####### Activity functions: Human movement ####
movementTimeStep = function(HPop, mobmat,day,control_df,prob_move_per_TS){
  movement_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  daynew = day
  while (length(which(mobmat$date == daynew)) == 0){
    daynew = daynew - days(1)
  }

  daily_move = as.data.frame(subset(mobmat,date == daynew))
  
  daily_move_mat = daily_move[,is.element(names(daily_move),c("fr_ID","to_ID","movenew"))]
  daily_move_mat = as.matrix(daily_move_mat)
  col1 = which(colnames(daily_move_mat) == "fr_ID")
  col2=which(colnames(daily_move_mat) == "to_ID")
  
  colmove = which(colnames(daily_move_mat) == "movenew")
  movement_matrix[daily_move_mat[,c(col1,col2)]] = daily_move_mat[,colmove]
  if (length(which(rowSums(movement_matrix)>1)) > 0){
   # print("Warning: row sums > 1 in movement matrix. Correcting...")
    correctingrows = which(rowSums(movement_matrix)>1)
    for (i in correctingrows){
    movement_matrix[i,] = movement_matrix[i,] /sum(movement_matrix[i,] )
    }
  }
  if (prob_move_per_TS > 0){
    movement_matrix = movement_matrix*prob_move_per_TS
  }
  
  diag(movement_matrix)=NA
  diag(movement_matrix) = 1-rowSums(movement_matrix,na.rm=T)
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
  HPop$nInf = round(HPop$nInf)
  
  # stochastic version
    z <- rbinom(n=NPat^2,size = rep(HPop$nInf,each=NPat),prob = t(movement_matrix)[])
  moved_matrix = t(matrix(z,NPat,NPat,dimnames=list(patIDs,patIDs)))
  for (i in 1:dim(moved_matrix)[1]){
     if (sum(moved_matrix[i,]) > 0){
     moved_matrix[i,] = moved_matrix[i,]/sum(moved_matrix[i,]) * HPop$nInf[i]
     }
  }
  diag(moved_matrix)=0
  #print(sum(moved_matrix))
  HPop$nInfMovedToday = colSums(moved_matrix)

  HPop$nInf = HPop$nInf - rowSums(moved_matrix) + colSums(moved_matrix)
  
  for (i in 1:length(HPop$nImportedTodayIntl)){
    HPop$nImportedTodayIntl[i] = sum(moved_matrix[which(HPop$country != HPop$country[i]),i])
    HPop$nExportedTodayIntl[i] = sum(moved_matrix[i,which(HPop$country != HPop$country[i])])
  }
  HPop$nImportedToday = colSums(moved_matrix)
  HPop$nExportedToday = rowSums(moved_matrix)
  #quick fix
  for (i in 1:length(HPop$nInf)){
    if (HPop$nInf[i] > HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec[i]){
     HPop$nRec[i] = max(0, HPop$nTotal[i] - HPop$nExp[i]- HPop$nInf[i])
      
      HPop$nInf[i] = max(0,HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec[i])
    }
    if (HPop$nInf[i] <0 ){
      
      HPop$nInf[i] = 0
    }
  }

  #(max((HPop$nInf + HPop$nRec + HPop$nExp)/HPop$nTotal))
  #print(paste0("Number of infected people moving: ",sum(abs(HPop$nInfMovedToday))/2))
  list(HPop=HPop,moved_matrix=moved_matrix)
}


###### Response functions: Control
#relative_movement is the proportion of original movement out/in that we want to keep -- ie. .1 = 10% of original movement rate
stopMovement = function(HPop,mobmat,current_date){

  stopping = which(HPop$controlled > 0)
  if (length(stopping) > 0){
    #print(paste("Stopping movement in ",length(HPop$names[stopping]), " of ", length(HPop$nInf)))
    for (ctrl_pat in stopping){
      control_patches = HPop$IDs[ctrl_pat]
      mobmat[control_patches,] = mobmat[control_patches,] * HPop$controlled[ctrl_pat]
      mobmat[,control_patches] = mobmat[,control_patches] * HPop$controlled[ctrl_pat]
    } 
    diag(mobmat)=NA
    if (sum(rowSums(mobmat,na.rm=T)>1) > 0 ){
      #print("Warning: row sums > 1 in movement matrix during stopMovement. Correcting...")
      correctingrows = which(rowSums(mobmat,na.rm=T)>1)
      for (i in correctingrows){
        print(sum(mobmat[i,]))
        
        mobmat[i,] = mobmat[i,] /sum(mobmat[i,] )
      }
      }
    diag(mobmat) = 1-rowSums(mobmat,na.rm=T)
  }
  mobmat
}


reduceContacts = function(HPop,contact_reduction_df,current_date){
  daynew = current_date
  while (length(which(contact_reduction_df$date == daynew)) == 0){
    daynew = daynew - days(1)
  }
  
  HPop$relativeContactRate = rep(1,length(HPop$names))
    contact_reduction_df_sub = subset(contact_reduction_df,date == daynew)
    if (dim(contact_reduction_df_sub)[1] > 0){
      #print(paste("Reducing contacts in ",length(contact_reduction_df_sub$from), " of ", length(HPop$nInf)))
      for (i in 1:dim(contact_reduction_df_sub)[1]){
        HPop$relativeContactRate[which(HPop$names == contact_reduction_df_sub$from[i])] = contact_reduction_df_sub$relative_contact[i]
        
      }
    }
  
  HPop
}



###### Master function  ####
runSim_new = function(HPop,pat_info,movement_reduction_dfnold, contact_reduction_dfnold,movement_reduction_dfld, contact_reduction_dfld,mobmat,day_list,recrate_values,exposerate_df,exposepd,exposed_pop_inf_prop = 0,TSinday = 1,prob_move_per_TS=0) {
  
  
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
      exposerate_df$exposerate = exposerate_df$exposerate/TSinday
    }
  }
  total_movement = matrix(0,length(HPop$nInf),length(HPop$nInf))
  all_spread = matrix(0,length(day_list),length(HPop$nInf))
  all_spread_total = matrix(0,length(day_list),length(HPop$nInf))
  all_spread_today = matrix(0,length(day_list),length(HPop$nInf))
  colnames(all_spread) = HPop$names
  colnames(all_spread_total) = HPop$names
  colnames(all_spread_today) = HPop$names
  #print(all_dates)
  for (current_day in 1:length(day_list)){
    for (current_TS in 1:TSinday){
    #print()
    
    HPop = recoveryTimeStep(HPop,recrate_values,day_list[current_day])
    HPop = exposedtoinfTimeStep(HPop,1/exposepd)
    HPop = reduceContacts(HPop,contact_reduction_df,day_list[current_day])
    HPop = exposedTimeStep(HPop,exposerate_df, day_list[current_day], exposed_pop_inf_prop)
    outmove = movementTimeStep(HPop,mobmat,day_list[current_day],movement_reduction_df,prob_move_per_TS)
    HPop = outmove$HPop
    today_movement = outmove$moved_matrix
    total_movement = total_movement + today_movement
    
    }
    print(paste(day_list[current_day],"inf: ",sum(HPop$nInf)," exp:",sum(HPop$nExp), "rec: ",sum(HPop$nRec)))
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = day_list[current_day], inf = sum(HPop$nInf)))
    all_spread[current_day,] = HPop$nInf
    all_spread_total[current_day,] = HPop$nInf + HPop$nRec
    all_spread_today[current_day,] = HPop$nInfectedToday
    
  }
  all_spread_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_2= cbind(all_spread_2,all_spread)
  all_spread_total_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_total_2= cbind(all_spread_total_2,all_spread_total)
  all_spread_today_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_today_2= cbind(all_spread_today_2,all_spread_today)
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread_2,all_spread_today = all_spread_today_2,today_movement = today_movement,total_movement = total_movement,all_spread_total = all_spread_total_2)
}

