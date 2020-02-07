####Model running code for BEARmod v.0.5


source("bearmod/bearmod_fx.R")
#reading in data
movement_data = read.table("baidu/Baidu_IP_flow_201411_201505.txt",sep="\t",header=T)
cell_user_data = read.table("baidu/IPusers_201411-201511.txt",sep="\t",header=T)
pop_data = read.csv("baidu/2012SHP_2010census_2011_2016CDCPop_adm2.csv")
relative_move_data = read.csv("../Downloads/parameter_relative_movement.csv",stringsAsFactors=F)
relative_move_data$date = date(relative_move_data$date)

#patch names are retained to more easily link back to the original dataset
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
movement_data$Date = ymd(movement_data$date)


#Initial parameters
NPat = length(patNames)
patnInf = rep(0,NPat)
patnExp = c(rep(0,NPat) )

recrate = 1/6 #daily probability of recovery
exposerate = 2.68/6# R0 of 2.68, 5.8 days till seeking treatment # How many people a single person potentially infects per day -- can be calculated from R0 estimate if you divide R0 by infectious period
exposepd = 5 # incubation period


#start infection in Wuhan
patnInf[which(patNames == 42010000)] = 5



pat_locator = merge(pat_locator,pop_data[,c("SHP_CITY_CODE","TOTAL_POP2015")],by.x="patNames",by.y="SHP_CITY_CODE")

#recovery rate variable
recover_df = data.frame(date = seq(from=min(movement_data$Date),to=max(movement_data$Date),by="days"),recrate = recrate)
recover_df$recrate[which(recover_df$date > "2015-2-18")] = 1/3

### Running the model  ####



HPop = InitiatePop(pat_locator,patnInf,patnExp)
###dates of simulation

input_dates = seq(date("2015-1-5"),date("2015-3-25"),by="days") 



###### Master function  ####
HPop_update = runSim(HPop,pat_locator,relative_move_data,movement_data, input_dates,recover_df, exposerate,exposepd)


##Plotting the results
newHPop = HPop_update$HPop
epidemic_curve = HPop_update$epidemic_curve
all_spread = HPop_update$all_spread
plot(epidemic_curve$Date,epidemic_curve$inf)


# 
# ####Many iterations
# output_df=matrix(0,0,340)
# for (run in 1:1000){
#   print(run)
# HPop_update = runSim(HPop,pat_locator,movement_data, input_dates, recrate, exposerate,exposepd,relative_movement = 1)
# 
# 
# newHPop = HPop_update$HPop
# epidemic_curve = HPop_update$epidemic_curve
# all_spread = HPop_update$all_spread
# output_df = rbind(output_df,t(all_spread[,dim(all_spread)[2]]))
# }
# 
# output_df2 = output_df > 0
# colnames(output_df2) = newHPop$names
# 
# 
# chn_shp$prob_CNY = 0
# for (i in 1:dim(chn_shp)[1]){
#   if (length(which(colnames(output_df2) == chn_shp$ZONECODE[i]))>0){
#     chn_shp$prob_CNY[i] = colSums(output_df2)[which(colnames(output_df2) == chn_shp$ZONECODE[i])]/1000
#   }}
# ggplot() +
#   geom_sf(chn_shp, mapping = aes(fill = prob_CNY) ) + scale_fill_distiller(palette="YlOrRd",direction=1)
# ####Animation code
# chn_shp
# 
# ####Animation code
# 
# #Numbers of people infected per patch, with IDs
# outputdf = data.frame(ID = newHPop$names)
# output_data = cbind(outputdf,all_spread)
# library(reshape2)
# melt_data = melt(output_data,id.vars="ID")
# 
# library(ggplot2)
# library(ggmap)
# library(rgeos)
# library("plyr")
# library("ggplot2")
# library("maptools")
# library(raster)
# library(igraph)
# library(rgdal)
# library(MASS)
# library(fossil)
# library(McSpatial)
# library(geosphere)
# library(ggrepel)
# library(hexbin)
# library(gganimate)
# library(viridis)
# library(sf)
# chn_shp = read_sf(dsn="baidu/ChinaShapefile2012",layer="dishi")
# library(ggplot2)
# library(gganimate)
# library(ggmap)
# library(maps)
# library(gapminder)
# 
# 
# library(animation)
# 
# 
# dev.control('enable')
# 
# oopts = ani.options(interval = 0.3)
# 
# 
# ani.options(oopts)
# ani.record(reset = TRUE)
# i=1
# 
# plots=list()
# melt_data$variable = as.numeric(melt_data$variable)
# for (i in 1:(as.numeric(date("2015-1-30")-date("2014-12-01")))){
# melt_data2= subset(melt_data,as.numeric(variable) ==i)
# chn_shp2 = merge(chn_shp,melt_data2,by.x="ZONECODE",by.y="ID",all.x=T)
# chn_shp2$value= round(chn_shp2$value)
# plots[[i]] = ggplot() +
#   geom_sf(chn_shp2,mapping=aes(),fill="light grey") +
#   geom_sf(subset(chn_shp2,value>0), mapping = aes(fill = value) ) +
#   scale_fill_distiller(palette="YlOrRd",direction=1,limits = c(.1,max(melt_data$value)),trans="log10",name="# inf")+
#   ggtitle(i)
#   #scale_fill_viridis(direction=-1,option="A",trans="log10")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),                                                                            
#         panel.background = element_blank(), axis.line = element_line(colour = "white")) 
# print(i)
# 
# }
# trace.animate <- function(plotter) { for (i in 1:length(plotter)) {
#   print(plotter[[i]])
# }}
# saveGIF(trace.animate(plots),interval=.3,movie.name="cov_model.gif",ani.width=800,ani.height=800)
