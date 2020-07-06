library(janitor)
library(lubridate)
library(raster)
library(rgdal)
library(data.table)
library(tidyverse)
library(sf)
library(exactextractr)
nuts3 <- st_read("NUTS_RG_01M_2016_4326_LEVL_3.shp/NUTS_RG_01M_2016_4326_LEVL_3.shp",stringsAsFactors=F) # admin 3 units

nuts3$IDs = 1:dim(nuts3)[1]
#population raster from WorldPop.org
popdata = raster("ppp_2020_1km_Aggregated.tif")

nuts3$pop <-exactextractr::exact_extract(popdata, nuts3, fun="sum")
nuts3_popdf=st_drop_geometry(nuts3[,c("IDs","pop")])
#Google 2019 NUTS3 dataset; forms the baseline movement patterns for our model
movement_data = fread("agg_epi_mobility_corona",stringsAsFactors=F)


#Account for smartphone market share across countries when linking Google and Vodafone data
smartphone_penetration = read.csv("smartphone_penetration2018.csv")
smartphone_penetration$smartphone = as.numeric(gsub("[\\%,]", "", smartphone_penetration$Smartphone.penetration))/100


iso2 = read.csv("iso2_csv.csv",stringsAsFactors=F)
iso2$Code[which(iso2$Code=="GB")] = "UK"
smartphone_penetration$ISO = ""
for (i in 1:dim(smartphone_penetration)[1]){
  if (length(which(iso2$Name == smartphone_penetration$Country[i]))){
  smartphone_penetration$ISO[i] = iso2$Code[which(iso2$Name == smartphone_penetration$Country[i])]
  }}
#this dataset has some islands and other areas that aren't a part of the mainland country, so we'll remove and just use the country-level estimate
smartphone_penetration=smartphone_penetration[-c(11,15,25,29,33,35),]
nuts3$penetration = 0
for (i in 1:dim(nuts3)[1]){
  if (length(which(smartphone_penetration$ISO == nuts3$CNTR_CODE[i]))){
    nuts3$penetration[i] = smartphone_penetration$smartphone[which(smartphone_penetration$ISO == nuts3$CNTR_CODE[i])]
  }}
nuts3$penetration[which(nuts3$penetration == 0)] = mean(nuts3$penetration[-which(nuts3$penetration ==0 )])
movement_data$to_name = substr(movement_data$V4,start=4,stop=8)
movement_data$fr_name = substr(movement_data$V3,start=4,stop=8)


movement_data$to_ID = 0
movement_data$fr_ID =0
movement_data$fr_penetration =0
movement_data$to_ID = apply(movement_data,1,FUN=function(x){nuts3$IDs[which(nuts3$NUTS_ID == x[6])]})
movement_data$fr_ID = apply(movement_data,1,FUN=function(x){nuts3$IDs[which(nuts3$NUTS_ID == x[7])]})

movement_data$fr_penetration = apply(movement_data,1,FUN=function(x){nuts3$penetration[which(nuts3$NUTS_ID == x[7])]})
#the movement value itself is the log of a relative flow value. so we'll exponentiate to get the original relative population flow value
movement_data$move = exp(movement_data$V5)

movement_data$pop_fr = 0
movement_data$pop_fr = apply(movement_data,1,FUN=function(x){nuts3$pop[which(nuts3$NUTS_ID == x[7])]})
#we adjust population flow values usign population in NUTS3area, further adjusted for smartphone penetration in that area.
movement_data$glh = log(movement_data$move/(movement_data$pop_fr * movement_data$fr_penetration))
#this linear model includes only "glh" as a predictor for the probability of transition in the vodafone data
movement_model = readRDS("glm_model_log.gz")

#predict the probability of transition using the adjusted google flow
movement_data$movenew = exp(predict(movement_model,re.form=NA,newdata=movement_data))

patNames = nuts3$NUTS_ID
patIDs = nuts3$IDs
pat_locator = data.frame(patNames,patIDs,nuts3$pop)
names(pat_locator)[which(names(pat_locator)== "nuts3.pop")] = "pop"
pat_locator$pop[which(is.na(pat_locator$pop))] = 1


