
library(lubridate)

# ## 2013 - 2014 data
# movement_data2 = read.csv("baidu/7 feb/Baidu_LBS_flow_201312-201404.csv")
# cell_user_from_data = read.csv("baidu/7 feb/LBSusers_from_201312-201404.csv")
# cell_user_from_data$date = date(cell_user_from_data$date) + days(1)
# cell_user_to_data = read.csv("baidu/7 feb/LBSusers_to_201312-201404.csv")
# cell_user_to_data$date = date(cell_user_to_data$date) + days(1)
pop_data = read.csv("baidu/2012SHP_2010census_2011_2016CDCPop_adm2.csv")

relative_move_data = read.csv("baidu/7 feb/parameter_relative_movement_2014 8 feb.csv",stringsAsFactors=F)
relative_move_data$date = date(relative_move_data$date) + years(1) + days(19)
#### 2015 data

movement_data = read.table("baidu/Baidu_IP_flow_201411_201505.txt",sep="\t",header=T)
cell_user_data = read.table("baidu/IPusers_201411-201511.txt",sep="\t",header=T)
pop_data = read.csv("baidu/2012SHP_2010census_2011_2016CDCPop_adm2.csv")
###

patNames = unique(movement_data$to)[order(unique(movement_data$to))]  
patIDs = 1:length(patNames)
pat_locator = data.frame(patNames,patIDs)

#add patch names to the movement data
movement_data = merge(movement_data,pat_locator,by.x = "from",by.y = "patNames")
names(movement_data)[which(names(movement_data) == "patIDs")] = "fr_pat"
movement_data = merge(movement_data,pat_locator,by.x = "to",by.y = "patNames")
names(movement_data)[which(names(movement_data) == "patIDs")] = "to_pat"

#add total numbers of users to movement data
movement_data = merge(movement_data,cell_user_data,by.x = c("from","date"),by.y=c("SHP_CITY_CODE","date"))
names(movement_data)[which(names(movement_data) == "users")] = "fr_users"

movement_data = merge(movement_data,cell_user_data,by.x = c("to","date"),by.y=c("SHP_CITY_CODE","date"))
names(movement_data)[which(names(movement_data) == "users")] = "to_users"

#convert dates to format R can read
movement_data$date = ymd(movement_data$date)


names(movement_data)[which(names(movement_data) == "move")] = "movers"
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