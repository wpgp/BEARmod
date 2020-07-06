rm(list=ls())


#This runs the intermittent lockdown simulation runs
library(data.table) 
library(lubridate)
library(countrycode)
library(parallel)
library(countrycode)

source("eu_model_init.R")

HPop = InitiatePop(pat_locator,patnInf,patnExp)
#we can stipulate if the lockdown effect is driven by data ("real"), or entirely simulated using parameters that determine relative contact when countries are under lockdown "simulated"
run_lockdown_type = "real"
run_lockdown_length = 4
run_lockdown_relative_contact = .35
run_no_lockdown_relative_contact = .8 #reduction in contact assumed if no lockdown occurred in "simulated" runs

countries = c("DE","UK","PL")

pivot_date = as.Date("2020-04-08")
end_lockdown_date = pivot_date + weeks(run_lockdown_length)

for (country_id in 1:length(countries)){
  country_left_out = countries[country_id]
glhflows_final = glhflows_all

if (run_lockdown_type == "simulated"){
  glhflows_final$r0 = mean(glhflows_final$r0)
  glhflows_lockdown$r0 = mean(glhflows_lockdown$r0)
  glhflows_lockdown_extra$r0 = mean(glhflows_lockdown_extra$r0)
  glhflows_nolockdown$r0 = mean(glhflows_nolockdown$r0)
}
while (max(glhflows_final$date) < end_lockdown_date){
 
  glhflows_add = glhflows_lockdown_extra
  if (run_lockdown_type == "simulated"){
    glhflows_add$relative_move = run_lockdown_relative_contact
    glhflows_add$relative_contact = glhflows_add$r0 * recprob *rescale_contact* run_lockdown_relative_contact
    if (!is.na(country_left_out)){

      glhflows_add$relative_move[glhflows_add$country == country_left_out] =run_no_lockdown_relative_contact
      glhflows_add$relative_contact[glhflows_add$country == country_left_out] = glhflows_add$r0[glhflows_add$country == country_left_out] * recprob *rescale_contact *run_no_lockdown_relative_contact
    }
  }
  
  if (run_lockdown_type == "real"){
    if (!is.na(country_left_out)){
      for (rowid in which(glhflows_add$country == country_left_out)){
      glhflows_add$relative_move[rowid] = glhflows_nolockdown$change[which(glhflows_nolockdown$NUTS3fr == glhflows_add$NUTS3fr[rowid])]
      glhflows_add$relative_contact[rowid] = glhflows_add$r0[rowid] * recprob *rescale_contact * glhflows_nolockdown$change[which(glhflows_nolockdown$NUTS3fr == glhflows_add$NUTS3fr[rowid])]
      }
    }
  }
  glhflows_add$date = max(glhflows_final$date) + days(7)
  #print(end_date - max(glhflows_final$date))
  glhflows_final = rbind(glhflows_final,glhflows_add)
}


while (max(glhflows_final$date) < end_date){
  
  glhflows_add = glhflows_nolockdown
  if (run_lockdown_type == "simulated"){
    glhflows_add$relative_move = run_no_lockdown_relative_contact
    glhflows_add$relative_contact = glhflows_add$r0 * recprob *rescale_contact * run_no_lockdown_relative_contact
  }
  glhflows_add$date = max(glhflows_final$date) + days(7)
  print(end_date - max(glhflows_final$date))
  glhflows_final = rbind(glhflows_final,glhflows_add)
}
exposerate = 1

glhflows_final$from = glhflows_final$NUTS3fr


HPop_new <- runSimbinder(glhflows_final)

HPop_new$epidemic_curve$lockdown_type = run_lockdown_type
HPop_new$epidemic_curve$sim_relative_contact = run_lockdown_relative_contact
HPop_new$epidemic_curve$lockdown_length = run_lockdown_length
HPop_new$epidemic_curve$country = countries[country_id]

HPop_new$all_spread$lockdown_type = run_lockdown_type
HPop_new$all_spread$sim_relative_contact = run_lockdown_relative_contact
HPop_new$all_spread$lockdown_length = run_lockdown_length
HPop_new$all_spread$country = countries[country_id]

all_movement = as.data.frame(HPop_new$total_movement)
all_movement$lockdown_type = run_lockdown_type
all_movement$sim_relative_contact = run_lockdown_relative_contact
all_movement$lockdown_length = run_lockdown_length
all_movement$country =countries[country_id]

all_spread_total = as.data.frame(HPop_new$all_spread_total)
all_spread_total$lockdown_type = run_lockdown_type
all_spread_total$sim_relative_contact = run_lockdown_relative_contact
all_spread_total$lockdown_length = run_lockdown_length
all_spread_total$country = countries[country_id]

basefilename = paste0("ldtype",run_lockdown_type,"_ldlength",run_lockdown_length,"_sim-contact",run_lockdown_relative_contact*10,"_ctry",countries[country_id])

write.csv(HPop_new$epidemic_curve,paste0("results/epicurve2_",basefilename,".csv"))
write.csv(HPop_new$all_spread,paste0("results/allspread2_",basefilename,".csv"))
write.csv(all_spread_total,paste0("results/all_spread_total2_",basefilename,".csv"))

}
