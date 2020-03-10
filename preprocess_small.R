
library(lubridate)

# ## 2013 - 2014 data
# movement_data2 = read.csv("baidu/7 feb/Baidu_LBS_flow_201312-201404.csv")
# cell_user_from_data = read.csv("baidu/7 feb/LBSusers_from_201312-201404.csv")
# cell_user_from_data$date = date(cell_user_from_data$date) + days(1)
# cell_user_to_data = read.csv("baidu/7 feb/LBSusers_to_201312-201404.csv")
# cell_user_to_data$date = date(cell_user_to_data$date) + days(1)

movement_data = read.table("testmove.csv",sep=",",header=T)

patNames = unique(movement_data$to)[order(unique(movement_data$to))]  
patIDs = 1:length(patNames)
pat_locator = data.frame(patNames,patIDs)

#convert dates to format R can read
movement_data$date = ymd("2020-05-01")

# 
# missing_dates = c(date("2014-1-17"), date("2014-2-2"),date("2014-2-18"),date("2014-2-20"),date("2014-3-1"),date("2014-3-2"))
# for (dates in 1:length(missing_dates)){
#   replaceday = subset(movement_data,Date == missing_dates[dates] - days(1))
#   replaceday$Date = replaceday$Date + days(1)
#   movement_data = rbind(movement_data,replaceday)
# }

recrate = 1/6 #daily probability of recovery
exposerate = 2.68/6 # R0 of 2.68, 5.8 days till seeking treatment # How many people a single person potentially infects per day -- can be calculated from R0 estimate if you divide R0 by infectious period
exposepd = 3 # incubation period