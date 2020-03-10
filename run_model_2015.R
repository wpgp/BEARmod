####Model running code for BEARmod v.0.6
rm(list=ls())
library(data.table) # fread - fastly reading data
library(lubridate)

# setwd('//worldpop.files.soton.ac.uk/Worldpop/Projects/WP519091_Seasonality')
# setwd('D:/OneDrive - University of Southampton/Wuhan Coronavirus R0/Spread risk')
#setwd('C:/Users/sl4m18/OneDrive - University of Southampton/Wuhan Coronavirus R0/Spread risk')


source("bearmod/BEARmod_development/bearmod_fx_dev.R")
# source("bearmod/bearmod_fx.R")
source("bearmod/BEARmod_development/preprocess_data_2015_dev.R")
#Initial parameters
NPat = length(patNames)
patnInf = rep(0,NPat)
patnExp = c(rep(0,NPat) )



#start infection in Wuhan
patnInf[which(patNames == 42010000)] = 10


# pop2014 or pop2015
pat_locator = merge(pat_locator,pop_data[,c("SHP_CITY_CODE","TOTAL_POP2015")],by.x="patNames",by.y="SHP_CITY_CODE")
names(pat_locator)[which(names(pat_locator) == "TOTAL_POP2015")] = "pop"

#recovery rate variable
recover_df = data.frame(date = seq(from=min(movement_data$date),to=max(movement_data$date),by="days"),recrate = recrate)
 recover_df$recrate[which(recover_df$date > "2015-2-18")] = 1/3
 
 
#### Running the model  ####



HPop = InitiatePop(pat_locator,patnInf,patnExp)
###dates of simulation

# input_dates = seq(date("2015-1-5"),date("2015-3-25"),by="days") 
# input_dates = seq(date("2015-1-2"),date("2015-3-3"),by="days") # from 2020-12-08 to 2 wks after LNY's day 
# input_dates = seq(date("2015-1-2"),date("2015-3-17"),by="days") # coresponding to the period from 2020-12-08 to 4 wks after LNY's day 
#day_list = seq(date("2013-12-02"),date("2014-2-13"),by="days")

 input_dates = seq(date("2014-12-26"),date("2015-2-1"),by="days")
 # input_dates = seq(date("2013-12-02"),date("2014-2-13"),by="days") # coresponding to the period from 2020-12-08 to 2 wks after LNY's day 
# input_dates = seq(date("2013-12-02"),date("2014-2-27"),by="days") # coresponding to the period from 2020-12-08 to 4 wks after LNY's day
results = list()

for (run in 1:500){
  
  HPop_update = runSim(HPop,pat_locator,relative_move_data,movement_data, input_dates,recover_df, exposerate,exposepd,exposed_pop_inf_prop = .25, TSinday = 1)
  print(paste0("Run # ",run))
  results[[run]] = HPop_update$all_spread
}
save(results,file="results.RData")
# 
# ###### Master function  ####
# HPop_update = runSim(HPop,pat_locator,relative_move_data,movement_data, input_dates,recover_df, exposerate,exposepd)
# 
# 
# ##Plotting the results
# newHPop = HPop_update$HPop
# epidemic_curve = HPop_update$epidemic_curve
# all_spread = HPop_update$all_spread
# epidemic_curve$Date_2020 <- seq(as.Date('2020-12-08'),(as.Date('2020-12-08') + (max(input_dates) - min(input_dates))), by='days') 
# epidemic_curve$Date_2020_report <- epidemic_curve$Date_2020 + exposepd + 4 # 4 days - diagnosis, report delay
# epidemic_curve$inf_accu <- epidemic_curve$inf 
# for(i in 2:nrow(epidemic_curve)){
#   epidemic_curve$inf_accu[i] <- epidemic_curve$inf_accu[i-1] + epidemic_curve$inf[i]   
# }
# 
# plot(epidemic_curve$Date_2020,epidemic_curve$inf)
# plot(epidemic_curve$Date_2020_report,epidemic_curve$inf)
# plot(epidemic_curve$Date_2020_report,epidemic_curve$inf_accu)
# 
# 
# # 
# # ####Many iterations
# # output_df=matrix(0,0,340)
# # for (run in 1:1000){
# #   print(run)
# # HPop_update = runSim(HPop,pat_locator,movement_data, input_dates, recrate, exposerate,exposepd,relative_movement = 1)
# # 
# # 
# # newHPop = HPop_update$HPop
# # epidemic_curve = HPop_update$epidemic_curve
# # all_spread = HPop_update$all_spread
# # output_df = rbind(output_df,t(all_spread[,dim(all_spread)[2]]))
# # }
# # 
# # output_df2 = output_df > 0
# # colnames(output_df2) = newHPop$names
# # 
# # 
# # chn_shp$prob_CNY = 0
# # for (i in 1:dim(chn_shp)[1]){
# #   if (length(which(colnames(output_df2) == chn_shp$ZONECODE[i]))>0){
# #     chn_shp$prob_CNY[i] = colSums(output_df2)[which(colnames(output_df2) == chn_shp$ZONECODE[i])]/1000
# #   }}
# # ggplot() +
# #   geom_sf(chn_shp, mapping = aes(fill = prob_CNY) ) + scale_fill_distiller(palette="YlOrRd",direction=1)
# # ####Animation code
# # chn_shp
# # 
# # ####Animation code
# # 
# # #Numbers of people infected per patch, with IDs
# # outputdf = data.frame(ID = newHPop$names)
# # output_data = cbind(outputdf,all_spread)
# # library(reshape2)
# # melt_data = melt(output_data,id.vars="ID")
# # 
# # library(ggplot2)
# # library(ggmap)
# # library(rgeos)
# # library("plyr")
# # library("ggplot2")
# # library("maptools")
# # library(raster)
# # library(igraph)
# # library(rgdal)
# # library(MASS)
# # library(fossil)
# # library(McSpatial)
# # library(geosphere)
# # library(ggrepel)
# # library(hexbin)
# # library(gganimate)
# # library(viridis)
# # library(sf)
# # chn_shp = read_sf(dsn="baidu/ChinaShapefile2012",layer="dishi")
# # library(ggplot2)
# # library(gganimate)
# # library(ggmap)
# # library(maps)
# # library(gapminder)
# # 
# # 
# # library(animation)
# # 
# # 
# # dev.control('enable')
# # 
# # oopts = ani.options(interval = 0.3)
# # 
# # 
# # ani.options(oopts)
# # ani.record(reset = TRUE)
# # i=1
# # 
# # plots=list()
# # melt_data$variable = as.numeric(melt_data$variable)
# # for (i in 1:(as.numeric(date("2015-1-30")-date("2014-12-01")))){
# # melt_data2= subset(melt_data,as.numeric(variable) ==i)
# # chn_shp2 = merge(chn_shp,melt_data2,by.x="ZONECODE",by.y="ID",all.x=T)
# # chn_shp2$value= round(chn_shp2$value)
# # plots[[i]] = ggplot() +
# #   geom_sf(chn_shp2,mapping=aes(),fill="light grey") +
# #   geom_sf(subset(chn_shp2,value>0), mapping = aes(fill = value) ) +
# #   scale_fill_distiller(palette="YlOrRd",direction=1,limits = c(.1,max(melt_data$value)),trans="log10",name="# inf")+
# #   ggtitle(i)
# #   #scale_fill_viridis(direction=-1,option="A",trans="log10")+
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),                                                                            
# #         panel.background = element_blank(), axis.line = element_line(colour = "white")) 
# # print(i)
# # 
# # }
# # trace.animate <- function(plotter) { for (i in 1:length(plotter)) {
# #   print(plotter[[i]])
# # }}
# # saveGIF(trace.animate(plots),interval=.3,movie.name="cov_model.gif",ani.width=800,ani.height=800)
