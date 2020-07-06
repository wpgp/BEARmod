

runSimbinder = function(glhflows_input){
  #3 times steps per day (8 hours), 40% chance of moving per time step
  HPop_out =  runSim(HPop,pat_locator,glhflows_input,glhflows_input,movement_data, input_dates,recover_df, exposerate,exposepd,exposed_pop_inf_prop = 0, TSinday = 3,prob_move_per_TS=.4)
  HPop_out
}

source("bearmod_fx_dev_EU_version.R")
source("preprocess_data_EU.R")
#Initial parameters
NPat = length(patNames)
patnInf = rep(0,NPat)
patnExp = c(rep(0,NPat) )

pat_locator$country = substr(pat_locator$patNames,1,2)

r0vals = read.csv("pop_data_NPI_startDate_v2.csv")

recprob = 1/4
recrate = 1-exp(-1/4) #daily probability of recovery
exposerate = 1 # This is set to 1 because we will set exposure rates on a NUTS3-level basis (ie. accounting for country-level heterogeneity in R0)

pat_locator$ISO3 = countrycode(pat_locator$country,origin="eurostat",destination="iso3c")
pat_locator$r0 = 0
#define country-level R0
for (i in 1:dim(pat_locator)[1]){
  pat_locator$r0[i] = subset(r0vals,iso3c ==pat_locator$ISO3[i])$R0.ML
}

#read in cases on March 20 across Europe
patinfdata = read.csv("nuts3_cases_320.csv",stringsAsFactors=F)

#define cases based on cases on march 20 across all NUTS3 areas
for (i in 1:length(patNames)){
  patnInf[which(patinfdata$NUTS_ID == patNames[i])] = patinfdata$nuts3_cases[i]
}
library(viridis)


#Repeat mobility data to add 6 months to the end of the  Google NUTS3 dataset
movement_data$weekyr = paste("1",movement_data$V2+1,movement_data$V1,sep="-")
movement_data=subset(movement_data,V2<52)
movement_data$dateorig = as.Date(movement_data$weekyr,format="%w-%W-%Y")
movement_data$date = movement_data$dateorig  + months(2)

movement_data_add = subset(movement_data,date>"2019-03-08")
movement_data_add$date = movement_data_add$date + months(3)
movement_data = rbind(movement_data,movement_data_add)
movement_data_add$date = movement_data_add$date + months(3)
movement_data = rbind(movement_data,movement_data_add)
#Because the Google NUTS3 dataset is from 2019, we will add one year to simulate 2020
movement_data$date = movement_data$date + years(1)

start_date = as.Date("2020-03-20")
movement_data = subset(movement_data,date>(start_date - weeks(1)))
end_date = as.Date("2020-03-20") + weeks(30)
movement_data = subset(movement_data,date < end_date)
input_dates=seq(start_date,end_date,by="days")
recover_df = data.frame(date = input_dates,recrate = recprob)

#This dataset records the relative movement in a given week, compared to a baseline, for 2020 (Recording reductions in movement due to COVID-19)
glhflows = read.csv("glh flow data_NUTS3_EU.csv",stringsAsFactors=F)
glhflows = glhflows[c("NUTS3fr","t1","change")]
glhflows_fillin = expand.grid(NUTS3fr = pat_locator$patNames[!is.element(pat_locator$patNames,glhflows$NUTS3fr)], t1 = unique(glhflows$t1) )
glhflows_fillin$change = 0 
#for areas that did not have a reduction recorded, we used the mean value observed across all NUTS3 areas
#note: t1 represents the first day of a given week -- ie. "01/03/2020" means the movement observed from March 1 to March 7
for (date in levels(glhflows_fillin$t1)){
  print(date)
  glhflows_fillin$change[glhflows_fillin$t1 == date] = mean(glhflows$change[glhflows$t1 == date])
}

glhflows_all = rbind(glhflows,glhflows_fillin)
glhflows_all = merge(glhflows_all,pat_locator[c("patNames","r0","country")],by.x="NUTS3fr",by.y="patNames")
glhflows_all$change[is.na(glhflows_all$change)] = mean(subset(glhflows_all,t1 == "08/03/2020")$r0,na.rm=T)
#new NUTS3-specific contact rate calculated by also accounting for locak R value
glhflows_all$relative_contact = glhflows_all$change * glhflows_all$r0 * recrate 
glhflows_all$relative_move = glhflows_all$change
glhflows_all$date = dmy(glhflows_all$t1)
#NUTS3 areas with no lockdown inherited the movement patterns from March 1
glhflows_nolockdown = subset(glhflows_all,t1=="01/03/2020")



glhflows_lockdown = subset(glhflows_all,t1=="22/03/2020")
#We extended the lockdown effect one further week into April, to account for additional community uptake of lockdown
glhflows_lockdown_extra = subset(glhflows_all,t1=="22/03/2020")
glhflows_lockdown_extra$change =glhflows_lockdown_extra$change *.75
#Both contact rates and movement reduced based on reductions in movement observed in the Google COVID-19 dataset
glhflows_lockdown_extra$relative_contact = glhflows_lockdown_extra$change * glhflows_lockdown_extra$r0 * recrate 
glhflows_lockdown_extra$relative_move = glhflows_lockdown_extra$change 

