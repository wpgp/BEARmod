library(janitor)
library(lubridate)
library(raster)

ita_adm3 = readOGR(dsn="gadm36_ITA_shp",layer="gadm36_ITA_3")
ita_adm3$IDs = 1:dim(ita_adm3)[1]
ita_adm3$NAME_3_lower = tolower(ita_adm3$NAME_3)

#Movement matrix
movement_data = read.csv("odvals.csv",stringsAsFactors=F)


movement_data$from = tolower(movement_data$origin)
movement_data$to = tolower(movement_data$commune)
movement_data$date = as.Date(movement_data$date)

cell_user_data$commune_lower = movement_data$from

#User number matrix
cell_user_data = read.csv("allvals.csv",stringsAsFactors=F)


cell_user_data$commune_lower = tolower(cell_user_data$commune)
cell_user_data$date = as.Date(cell_user_data$date)
cell_user_data = cell_user_data[,c("date","commune_lower","distinct","all")]

###
patNames = ita_adm3$NAME_3_lower
patIDs = ita_adm3$IDs
pat_locator = data.frame(patNames,patIDs)

#add patch names to the movement data
movement_data = merge(movement_data,pat_locator,by.x = "from",by.y = "patNames")
names(movement_data)[which(names(movement_data) == "patIDs")] = "fr_pat"
movement_data = merge(movement_data,pat_locator,by.x = "to",by.y = "patNames")
names(movement_data)[which(names(movement_data) == "patIDs")] = "to_pat"

#add total numbers of users to movement data
movement_data = merge(movement_data,cell_user_data,by.x = c("from","date"),by.y=c("commune_lower","date"))
names(movement_data)[which(names(movement_data) == "distinct")] = "fr_users"

movement_data = merge(movement_data,cell_user_data,by.x = c("to","date"),by.y=c("commune_lower","date"))
names(movement_data)[which(names(movement_data) == "distinct")] = "to_users"


names(movement_data)[which(names(movement_data) == "movs")] = "movers"
# 
recrate = 1/6 #daily probability of recovery
exposerate = 2.68/6 # R0 of 2.68, 5.8 days till seeking treatment # How many people a single person potentially infects per day -- can be calculated from R0 estimate if you divide R0 by infectious period
exposepd = 3 # incubation period