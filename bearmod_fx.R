##### BEARmod v.0.92
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
InitiatePop = function(pat_locator,initialInf,initialExp,asymp_frac = 0){
  NPat = dim(pat_locator)[1]
  list(
    nInitialInf = initialInf,
    nInitialExp = initialExp,
    nInf_symp = initialInf * (1 - asymp_frac),
    nInf_asymp = initialInf * asymp_frac,
    nQ_Fatal = rep(0,NPat),
    nQ_notFatal = rep(0,NPat),
    nExp = initialExp,
    nRec_symp = rep(0,NPat),
    nRec_asymp = rep(0,NPat),
    nDeath = rep(0,NPat),
    nTotal = pat_locator$pop,
    names = pat_locator$patNames,
    IDs = pat_locator$patIDs,
    relativeInf = rep(1,NPat),
    nRecoveredToday = rep(0,NPat),
    nInfectedToday_symp = rep(0,NPat),
    nInfectedToday_asymp = rep(0,NPat),
    nExposedToday = rep(0,NPat),
    nInfMovedToday = rep(0,NPat),
    controlled = rep(0,NPat)
  )
}

##### Epidemic functions: exposure, infectivity, recovery  ####
recoveryTimeStep = function(HPop, recrate_values,current_day){
  symp_recrate = subset(recrate_values,date == current_day)$recrate_symp
  asymp_recrate = subset(recrate_values,date == current_day)$recrate_asymp
  HPop$nInf_asymp = round(HPop$nInf_asymp)
  HPop$nQ_notFatal = round(HPop$nQ_notFatal)
 #print(recrate)#print(paste0("Day ",current_day, " recovery rate: ", recrate))
  for (i in 1:length(HPop$nInf_asymp)){
    
    recQ_symp = sum(rbinom(HPop$nQ_notFatal[i],1,symp_recrate))
    recInf_asymp = sum(rbinom(HPop$nInf_asymp[i],1,asymp_recrate))
    HPop$nRecoveredToday[i]= recQ_symp +recInf_asymp
    HPop$nInf_asymp[i] = HPop$nInf_asymp[i] - recInf_asymp
    HPop$nQ_notFatal[i] = HPop$nQ_notFatal[i] - recQ_symp
    
    HPop$nRec_symp[i] = HPop$nRec_symp[i] + recQ_symp
    HPop$nRec_asymp[i] = HPop$nRec_asymp[i] +recInf_asymp
    
  }
  #print(paste0("Number of people recovering: ",sum(HPop$nRecoveredToday)))
  HPop
}

quarantineTimeStep = function(HPop, quarantine_df,current_day,probability_death){
  daynew = current_day
  
  while (length(which(recrate_values$date == daynew)) == 0){
    daynew = daynew - days(1)
  }
  
  daily_move = as.data.frame(subset(mobmat,date == daynew))
  
  quarantine_prob = subset(quarantine_df,date == daynew)$quarantine_prob
  
  
  nQ_Fatal = rep(0,NPat),
  nQ_notFatal = rep(0,NPat),
  HPop$nInf_symp = round(HPop$nInf_symp)
  for (i in 1:length(HPop$nInf_symp)){
    
    inf_progressed_to_quarantine = sum(rbinom(HPop$nInf_symp[i],1,quarantine_prob))
    quarantine_to_die = sum(rbinom(inf_progressed_to_quarantine,1,probability_death))
    quarantine_not_to_die = inf_progressed_to_quarantine - quarantine_to_die
    
    
    HPop$nInf_asymp[i] = HPop$nInf_asymp[i] - recInf_asymp
    HPop$nQ_notFatal[i] = HPop$nQ_notFatal[i] - recQ_symp
    
    HPop$nRec_symp[i] = HPop$nRec_symp[i] + recQ_symp
    HPop$nRec_asymp[i] = HPop$nRec_asymp[i] +recInf_asymp
    
  }
  #print(paste0("Number of people recovering: ",sum(HPop$nRecoveredToday)))
  HPop
}

exposedtoinfTimeStep = function(HPop, exp_to_infrate,prop_asymptomatic = 0){
  #(exp_to_infrate)
  HPop$nExp = round(HPop$nExp)
  for (i in 1:length(HPop$nInf_symp)){
    total_infected_today = sum(rbinom(HPop$nExp[i],1,exp_to_infrate))
    HPop$nInfectedToday_asymp[i]= HPop$nInfectedToday_asymp[i] * prop_asymptomatic
    
    HPop$nInfectedToday_symp[i]= HPop$nInfectedToday_symp[i] * (1-prop_asymptomatic)
   # }
    HPop$nExp[i] = HPop$nExp[i] - total_infected_today
    
  }
  HPop
}

exposedTimeStep = function(HPop, exposerate_df, current_day, exposed_pop_inf_prop){
  
  if (is.numeric(exposerate_df)){
    exposerate = exposerate_df
  }
  if (is.data.frame(exposerate_df)){
    exposerate = subset(exposerate_df, date == current_day)$exposerate
  }
  for (i in 1:length(HPop$nInf_symp)){ 
    infectious_pop = HPop$nInf_symp[i] +HPop$nInf_asymp[i]+ exposed_pop_inf_prop * HPop$nExp[i]
    infectious_pop = round(infectious_pop)
    HPop$nExposedToday[i]= sum(rpois(infectious_pop,exposerate)) * (1 - min(1, ( (HPop$nInf_symp[i] +HPop$nInf_asymp[i] + HPop$nExp[i] + HPop$nRec[i]) / HPop$nTotal[i]) ))
    
    
    if (HPop$nExp[i] + HPop$nExposedToday[i] < HPop$nTotal[i] - HPop$nInf_symp[i] -HPop$nInf_asymp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i]) {
      HPop$nExp[i] = HPop$nExp[i] + HPop$nExposedToday[i]
    
    } else {
      HPop$nExposedToday[i] = max(0,HPop$nTotal[i] - HPop$nInf_symp[i] -HPop$nInf_asymp[i] - HPop$nExp[i]- HPop$nRec_symp[i] - HPop$nRec_asymp[i])
      
      HPop$nExp[i]  = max(0,HPop$nTotal[i] -HPop$nInf_symp[i] -HPop$nInf_asymp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i])
      
    }
  }
  #print(paste0("Number of people newly exposed: ",sum(HPop$nExposedToday)))
  HPop
}


exposedTimeStep_timespent = function(HPop, exposerate_df, current_day, exposed_pop_inf_prop,ts_data){
  
  TS_matrix = matrix(0,NPat,NPat,dimnames=list(patIDs,patIDs))
  daily_move = subset(ts_data,date == current_day)
  daily_move = subset(daily_move,!is.na(fr_pat) & !is.na(to_pat) & !is.na(fr_users) & !is.na(movers))
  
  daily_move_mat = daily_move[,is.element(names(daily_move),c("fr_pat","to_pat","fr_users","movers"))]
  daily_move_mat = as.matrix(daily_move_mat)
  col1 = which(colnames(daily_move_mat) == "fr_pat")
  col2=which(colnames(daily_move_mat) == "to_pat")
  
  colmove = which(colnames(daily_move_mat) == "movers")
  colusers=which(colnames(daily_move_mat) == "fr_users")
  TS_matrix[daily_move_mat[,c(col1,col2)]] = daily_move_mat[,colmove]/daily_move_mat[,colusers]
  if (length(which(rowSums(TS_matrix)>1)) > 0){
    print("Warning: row sums > 1 in movement matrix. Correcting...")
    correctingrows = which(rowSums(TS_matrix)>1)
    for (i in correctingrows){
      TS_matrix[i,] = TS_matrix[i,] /sum(TS_matrix[i,] )
    }
  }
  for (i in 1:length(patIDs)){
    TS_matrix[i,i] = 1 - sum(TS_matrix[i,-i])
  }
  
  if (is.numeric(exposerate_df)){
    exposerate = exposerate_df
  }
  if (is.data.frame(exposerate_df)){
    exposerate = subset(exposerate_df, date == current_day)$exposerate
  }
  movement_adjusted_infectious_prop = rep(0,length(HPop$nInf))
  for (i in 1:length(HPop$nInf)){
    movement_adjusted_infectious_prop[i] = sum(((HPop$nInf * TS_matrix[,i]) + exposed_pop_inf_prop * sum(( HPop$nExp * TS_matrix[,i])))) / sum(HPop$nTotal * TS_matrix[,i])
  }
  susceptible_vec = HPop$nTotal - HPop$nInf - HPop$nExp - HPop$nRec_asymp- HPop$nRec_symp
  
  probability_infection = 1-exp(-exposerate * movement_adjusted_infectious_prop)
  for (i in 1:length(HPop$nInf)){
    susceptible_weighted_pop = round(susceptible_vec[i]*TS_matrix[i,])
    HPop$nExposedToday[i] = sum(rbinom(length(susceptible_weighted_pop),size = susceptible_weighted_pop,prob=probability_infection))
    if (HPop$nExp[i] + HPop$nExposedToday[i] < HPop$nTotal[i] - HPop$nInf[i] - HPop$nRec_asymp[i]- HPop$nRec_symp[i] ) {
      HPop$nExp[i] = HPop$nExp[i] + HPop$nExposedToday[i]
      
    } else {
      if (HPop$nExp[i]< 0){print(HPop$nExp[i])}
      
      HPop$nExposedToday[i] = HPop$nTotal[i] - HPop$nInf[i] - HPop$nExp[i]- HPop$nRec_symp[i]- HPop$nRec_asymp[i]
      HPop$nExp[i]  = HPop$nTotal[i] - HPop$nInf[i] - HPop$nRec_asymp[i]- HPop$nRec_symp[i]
      if (HPop$nExp[i]< 0){print(HPop$nExp[i])}
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
  HPop$nInf_symp = round(HPop$nInf_symp)
  
  HPop$nInf_asymp = round(HPop$nInf_asymp)
  # stochastic version
    moved_symp <- rbinom(n=NPat^2,size = rep(HPop$nInf_symp,each=NPat),prob = t(movement_matrix)[])
  moved_matrix_symp = t(matrix(moved_symp,NPat,NPat,dimnames=list(patIDs,patIDs)))
  moved_asymp <- rbinom(n=NPat^2,size = rep(HPop$nInf_asymp,each=NPat),prob = t(movement_matrix)[])
  moved_matrix_asymp = t(matrix(moved_asymp,NPat,NPat,dimnames=list(patIDs,patIDs)))
  
  for (i in 1:dim(moved_matrix_symp)[1]){
     if (sum(moved_matrix_symp[i,]) > 0){
     moved_matrix_symp[i,] = moved_matrix_symp[i,]/sum(moved_matrix_symp[i,]) * HPop$nInf_symp[i]
     }
  }
  for (i in 1:dim(moved_matrix_asymp)[1]){
    if (sum(moved_matrix_asymp[i,]) > 0){
      moved_matrix_asymp[i,] = moved_matrix_asymp[i,]/sum(moved_matrix_asymp[i,]) * HPop$nInf_asymp[i]
    }
  }
  #print(sum(moved_matrix))
  #print(sum(HPop$nInf))
  diag(moved_matrix_asymp)=0
  diag(moved_matrix_symp)=0
  HPop$nInfMovedToday = colSums(moved_matrix_asymp) + colSums(moved_matrix_asymp)

  HPop$nInf_symp = HPop$nInf_symp - rowSums(moved_matrix_symp) + colSums(moved_matrix_symp)
  HPop$nInf_asymp = HPop$nInf_asymp - rowSums(moved_matrix_asymp) + colSums(moved_matrix_asymp)
  #print(max((HPop$nInf + HPop$nRec + HPop$nExp)/HPop$nTotal))
  #quick fix
  for (i in 1:length(HPop$nInf_symp)){
    if (HPop$nInf_symp[i] > HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i] - HPop$nInf_asymp[i]){
     HPop$nRec_symp[i] = max(0, HPop$nTotal[i] - HPop$nExp[i]- HPop$nInf_symp[i] - HPop$nInf_asymp[i] - HPop$nRec_asymp[i])
      
      HPop$nInf_symp[i] = HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i] - HPop$nInf_asymp[i]
    }
    if (HPop$nInf_symp[i] <0 ){
      
      HPop$nInf_symp[i] = 0
    }
  }
  
  for (i in 1:length(HPop$nInf_asymp)){
    if (HPop$nInf_asymp[i] > HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i] - HPop$nInf_symp[i]){
      HPop$nRec_asymp[i] = max(0, HPop$nTotal[i] - HPop$nExp[i]- HPop$nInf_asymp[i] - HPop$nInf_symp[i]- HPop$nRec_symp[i])
      
      HPop$nInf_asymp[i] = HPop$nTotal[i] - HPop$nExp[i] - HPop$nRec_symp[i] - HPop$nRec_asymp[i] - HPop$nInf_symp[i]
    }
    if (HPop$nInf_asymp[i] <0 ){
      
      HPop$nInf_asymp[i] = 0
    }
  }
  
  #(max((HPop$nInf + HPop$nRec + HPop$nExp)/HPop$nTotal))
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
runSim = function(HPop,pat_info,control_info,mobmat,day_list,recrate_values,exposerate_df,exposepd,exposed_pop_inf_prop = 0,TSinday = 1,prob_move_per_TS=0,prop_asymptomatic = 0) {
  
  
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
  
  all_spread = matrix(0,length(day_list),length(HPop$nInf_symp))
  all_spread_today = matrix(0,length(day_list),length(HPop$nInf_symp))
  colnames(all_spread) = HPop$names
  #print(all_dates)
  for (current_day in 1:length(day_list)){
    for (current_TS in 1:TSinday){
    print(day_list[current_day])
    
    HPop = recoveryTimeStep(HPop,recrate_values,day_list[current_day])
    HPop = exposedtoinfTimeStep(HPop,1/exposepd,prop_asymptomatic)

    
    HPop = exposedTimeStep(HPop,exposerate_df, day_list[current_day], exposed_pop_inf_prop)
    
    HPop = movementTimeStep(HPop,mobmat,day_list[current_day],control_info,prob_move_per_TS)
    
    print(paste("inf: ",sum(HPop$nInf_symp)," exp:",sum(HPop$nExp), "rec: ",sum(HPop$nRec)))
    
    }
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = day_list[current_day], inf = sum(HPop$nInf)))
    all_spread[current_day,] = HPop$nInf_symp
    all_spread_today[current_day,] = HPop$nInfectedToday_symp
    
  }
  all_spread_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_2= cbind(all_spread_2,all_spread)
  all_spread_today_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_today_2= cbind(all_spread_today_2,all_spread_today)
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread_2,all_spread_today = all_spread_today_2)
}


runSim_timespent = function(HPop,pat_info,control_info,TS_data,day_list,recrate_values,exposerate_df,exposepd,exposed_pop_inf_prop = 0,TSinday = 1) {
  
  
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
  
  all_spread = matrix(0,length(day_list),length(HPop$nInf))
  colnames(all_spread) = HPop$names
  #print(all_dates)
  for (current_day in 1:length(day_list)){
    for (current_TS in 1:TSinday){
      print(day_list[current_day])
      HPop = recoveryTimeStep(HPop,recrate_values,day_list[current_day])
      HPop = exposedtoinfTimeStep(HPop,1/exposepd)
      HPop = exposedTimeStep_timespent(HPop,exposerate_df, day_list[current_day], exposed_pop_inf_prop,TS_data)
    }
    #save(HPop,file=paste(current_day,".RData"))
    epidemic_curve = rbind(epidemic_curve,data.frame(Date = day_list[current_day], inf = sum(HPop$nInf)))
    all_spread[current_day,] = HPop$nInf
    
  }
  all_spread_2 = data.frame(dates = day_list,runday = 1:length(day_list))
  all_spread_2= cbind(all_spread_2,all_spread)
  list(HPop = HPop,epidemic_curve = epidemic_curve,all_spread=all_spread_2)
}


