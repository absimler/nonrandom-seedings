### Nonrandom treatment locations and sagebrush seeding success##########
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######
#############################################################

### PART 1) Data formatting and organization  #############

library(betareg)
library(raster)
library(rgdal)
library(here)
library(sp)
library(rgeos)
library(spatialEco)
library(tidyverse)
library(lubridate)
library(rstanarm)
library(MatchIt)
library(tidybayes)
library(bayesplot)
library(cowplot)
library(brms)
library(modelr)
library(geosphere)
library(spdep)
library(loo)
library(patchwork)
library(cowplot)
library(egg)
library(dplyr)
library(tidyr)
library(stringr)
library(parallel)
library(rstan)
library(brms)

set.seed(1999)


## a) LOAD DATA: #############

data <- read.csv("formatteddata/data_feb9.csv") # treated locs
data2 <- read.csv("formatteddata/data2_feb9.csv") #untreated locs
data <- rbind(data, data2)



### b) Format data ####

### Reorganize structure of sagebrush cover variables so that we just retain the information about sage cover prefire (1yr), the year of the fire, 5 yrs postfire, and 10 yrs postfire
data$uniqueID <- c(1:40000)

sage <- gather(data, BITyear, sagecov, "X1985":"X2010")
sage <-  mutate(sage, BITyear=str_sub(BITyear, 2, 5))
sage$BITyear <- as.numeric(sage$BITyear)
sage$Fire_Year <- as.numeric(as.character(sage$Fire_Year))
sage$TSF <- sage$BITyear - sage$Fire_Year


weather <- gather(data, weatheryear, 
                  weatherval, "X1984ppt.pr":"X2015tmean.tmean")
weather <-  mutate(weather, year=str_sub(weatheryear, 2, 5))%>%
            mutate(weather, 
                   weathervar = str_sub(weatheryear, 6, 9))
weather$year <- as.numeric(weather$year)
weather$Fire_Year <- as.numeric(as.character(weather$Fire_Year))
weather$TSF <- weather$year - weather$Fire_Year

#Only retain the sage cover estimates from years needed (10 years pre-fire for assessment of pre-treatment trends, and a subset of postfire timepoints). Values over 100 represent various versions of NA (see RCMAP), so these have been eliminated from the pool of possible pixels to be selected.
sage <- subset(sage, (TSF==-10 | TSF==-9 | 
                        TSF==-8 | TSF==-7 |
                        TSF==-6 | TSF==-5 | 
                        TSF==-4 | TSF==-3 |
                        TSF==-2 | TSF==-1 | 
                        TSF==0 | TSF==1| TSF==9|
                        TSF==2 | TSF==3 | TSF==4 | 
                        TSF==6 | TSF==7 |
                        TSF==5 | TSF==8 | TSF==10 | 
                        TSF==12 |TSF==15) & sagecov<100) 

weather <- subset(weather, (TSF==-10 | TSF==-9 | 
                              TSF==-8 | TSF==-7 |
                              TSF==-6 | TSF==-5 | 
                              TSF==-4 | TSF==-3 |
                              TSF==-2 | TSF==-1 | 
                              TSF==0 | TSF==1| TSF==9|
                              TSF==2 | TSF==3 | TSF==4 | 
                              TSF==6 | TSF==7 |
                              TSF==5 | TSF==8 | 
                              TSF==10 | TSF==12 |TSF==15)) 
ppt <- subset(weather, weathervar=="ppt.")
temp <- subset(weather, weathervar=="tmea")

ppt <- ppt %>%
  dplyr::select(-year, -weatheryear) %>%
  spread(TSF, weatherval) 
temp <- temp %>%
  dplyr::select(-year, -weatheryear) %>%
  spread(TSF, weatherval) 

sage <- sage %>%
  dplyr::select(-BITyear) %>%
  spread(TSF, sagecov) 

sage <- sage %>% #rename variables for TSF
  rename(covTSF10pf="-10") %>%
  rename(covTSF9pf="-9") %>%
  rename(covTSF8pf="-8") %>%
  rename(covTSF7pf="-7") %>%
  rename(covTSF6pf="-6") %>%
  rename(covTSF5pf="-5") %>%
  rename(covTSF4pf="-4") %>%
  rename(covTSF3pf="-3") %>%
  rename(covTSF2pf="-2") %>%
  rename(covTSF1pf="-1") %>%
  rename(covTSF0 = "0") %>%
  rename(covTSF1 = "1") %>%
  rename(covTSF2 = "2") %>%
  rename(covTSF3 = "3") %>%
  rename(covTSF4 = "4") %>%
  rename(covTSF5 = "5") %>%
  rename(covTSF6 = "6") %>%
  rename(covTSF7 = "7") %>%
  rename(covTSF8 = "8") %>%
  rename(covTSF9 = "9") %>%
  rename(covTSF10 = "10") %>% 
  rename(elevation = Western_US_30m_DEM_extended_clip) %>%
  rename(heatload = Western_US_30m_Exp_HeatLo) %>%
  rename(distroad = disttoroad_westUS)

## Subset to putative "sagebrush habitat", which we are defining as places that had at least 1% sagebrush cover in one of the 5 years pre-fire:
sage <- subset(sage, covTSF1pf>0|covTSF2pf>0|covTSF3pf>0|
                   covTSF4pf>0|covTSF5pf>0)



### same formatting for the temperature and precip variables:
temp <- temp %>%
  rename(temp10pf="-10") %>%
  rename(temp9pf="-9") %>%
  rename(temp8pf="-8") %>%
  rename(temp7pf="-7") %>%
  rename(temp6pf="-6") %>%
  rename(temp5pf="-5") %>%
  rename(temp4pf="-4") %>%
  rename(temp3pf="-3") %>%
  rename(temp2pf="-2") %>%
  rename(temp1pf="-1") %>%
  rename(temp0 = "0") %>%
  rename(temp1 = "1") %>%
  rename(temp2 = "2") %>%
  rename(temp3 = "3") %>%
  rename(temp4 = "4") %>%
  rename(temp5 = "5") %>%
  rename(temp6 = "6") %>%
  rename(temp7 = "7") %>%
  rename(temp8 = "8") %>%
  rename(temp9 = "9") %>%
  rename(temp10 = "10")


ppt <- ppt %>%
  rename(ppt10pf="-10") %>%
  rename(ppt9pf="-9") %>%
  rename(ppt8pf="-8") %>%
  rename(ppt7pf="-7") %>%
  rename(ppt6pf="-6") %>%
  rename(ppt5pf="-5") %>%
  rename(ppt4pf="-4") %>%
  rename(ppt3pf="-3") %>%
  rename(ppt2pf="-2") %>%
  rename(ppt1pf="-1") %>%
  rename(ppt0 = "0") %>%
  rename(ppt1 = "1") %>%
  rename(ppt2 = "2") %>%
  rename(ppt3 = "3") %>%
  rename(ppt4 = "4") %>%
  rename(ppt5 = "5") %>%
  rename(ppt6 = "6") %>%
  rename(ppt7 = "7") %>%
  rename(ppt8 = "8") %>%
  rename(ppt9 = "9") %>%
  rename(ppt10 = "10")

ppt <- dplyr::select(ppt, uniqueID, ppt10pf:ppt10)
temp <- dplyr::select(temp, uniqueID, temp10pf:temp10)

## Subset to only the variables we'll need for this analysis:
pixels <- dplyr::select(sage, x, y, treatment, 
                        ID, Fire_Year, Fire_Name, 
                        Hectares, Shape_Area, Prj_Name, 
                        Trt_ID, Trt_Type, TrtMajCat,
                        Trt_Start, Trt_End, Artem_Sp_S, 
                        ARTR_Sp_Se, Grazing_Re, pt.ids,
                        Level3Ecoregions, 
                        elevation, clay, sand, heatload, 
                        distroad, 
                        febaprppt.pr:novjantmin.tmmn, 
                        uniqueID:covTSF10)

# The Ecoregion indicators DO NOT match up with numbering system from full dataset because our raster is clipped: Reassign numbers to names here to avoid confusion:
pixels$Level3Ecoregions[pixels$Level3Ecoregions==1] <- "Columbia Plateau"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==5] <- "Eastern Cascades\nFoothills"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==6] <- "NW Great Plains"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==8] <- "Middle Rockies"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==9] <- "Blue Mountains"

pixels$Level3Ecoregions[pixels$Level3Ecoregions==10] <- "Idaho Batholith"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==11] <- "Snake River\nPlain"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==12] <-                "Northern Basin\nand Range"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==13] <-                "Wyoming Basin"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==14] <-                "Central Basin\nand Range"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==15] <-               "High Plains"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==16] <-               "Wasatch and\nUinta Mountains"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==18] <-               "Sierra Nevada"
pixels$Level3Ecoregions[pixels$Level3Ecoregions==20] <-               "Southern Rockies"
               
pixels$Level3Ecoregions <- as.factor(pixels$Level3Ecoregions)

pixels <- dplyr::left_join(pixels, ppt, by="uniqueID")
pixels <- dplyr::left_join(pixels,temp, by="uniqueID")

#Subset treated areas to only seedings:
pixels.treat <- subset(pixels, treatment==1)
pixels.treat <- filter(pixels.treat, 
                       TrtMajCat=="Aerial Seeding" | 
                          Trt_Type=="Aerial Seeding: Fixed Wing" |
                          Trt_Type=="Aerial Seeding: Rotor Wing" |
                          Trt_Type=="Ground Seeding" |
                          Trt_Type=="ATV Broadcast" | 
                          Trt_Type == "Greenstrip: Drill Seeding" |
                          Trt_Type=="Ground Seeding: Drill" |
                          Trt_Type == "Ground Seeding: Imprinting"|
                          Trt_Type=="Broadcast"| 
                         Trt_Type=="Hand Broadcast"|
                          Trt_Type=="Seeding" | 
                         TrtMajCat=="Seeding")

pixels.untreat <- subset(pixels, treatment==0)


### c) Randomly select treated and untreated pixels: #######

## List of covariates that require values for analysis:
covars <- c('covTSF1pf',  "clay", "sand", 'febaprppt.pr',  
            'febaprtmean.tmean', 
            'novjanppt.pr', "novjantmean.tmean", 
            'Shape_Area', 'distroad', 'heatload',  
            "elevation", 
            'Level3Ecoregions', "Fire_Year",
            'covTSF0', "covTSF1", "covTSF2", 
            "covTSF3", "covTSF4", "covTSF5",
            "covTSF6", "covTSF7", "covTSF8", 
            "covTSF9", "covTSF10", "ppt1pf",
            'ppt0', "ppt1", "ppt2", "ppt3", 
            "ppt4", "ppt5",
            "ppt6", "ppt7", "ppt8", "ppt9", "ppt10", 
            "temp1pf", 'temp0', "temp1", "temp2", 
            "temp3", "temp4", "temp5",
            "temp6", "temp7", "temp8", "temp9", "temp10",
            "uniqueID")


# Remove observations without complete data for all covariates required.
pixels.treat <- pixels.treat %>%  
  dplyr::select(treatment, one_of(covars))%>%
  na.omit()

pixels.untreat <- pixels.untreat %>%  
  dplyr::select(treatment, one_of(covars))%>%
  na.omit()

# Select 10,000 treated and untreated pixels for equal sample sizes.
set.seed(2020201)
pixels.treat <- sample_n(pixels.treat, 10000)
pixels.untreat <- sample_n(pixels.untreat, 10000)

pixels.sub <- rbind(pixels.treat, pixels.untreat)

pixels.sub$Treatment_status <- ifelse(pixels.sub$treatment==1, "Treated", "Untreated")
pixels.sub$wintersprppt <- pixels.sub$febaprppt.pr + pixels.sub$novjanppt.pr


### d) Data summary statistics #######

## Correlation between covariates

corrmat <- cor(pixels.sub[c('covTSF1pf', 
                               'covTSF0', 'febaprppt.pr', "sand", 
                               "clay",  'febaprtmean.tmean',
                            'novjanppt.pr', "novjantmean.tmean", 
                               'Shape_Area',  'distroad',
                               "elevation", "heatload")])


## Table of summary statistics for covariates:

summ.stats <- as.data.frame(t(bind_rows(
  pixels.sub %>%
  group_by(treatment) %>%
  summarise_at(vars(covTSF1pf:covTSF10, -Level3Ecoregions),
               mean, na.rm = TRUE) %>%
  mutate(stat= ifelse(treatment==1, "treated mean", "untreated mean")),
  
  pixels.sub %>%
  group_by(treatment) %>%
  summarise_at(vars(covTSF1pf:covTSF10, -Level3Ecoregions
                    ),sd, na.rm = TRUE)%>%
  mutate(stat= ifelse(treatment==1, "treated sd", "untreated sd")),
  
  pixels.sub %>%
  group_by(treatment) %>%
  summarise_at(vars(covTSF1pf:covTSF10, 
                    -Level3Ecoregions),min, na.rm = TRUE) %>%
  mutate(stat= ifelse(treatment==1, "treated min", "untreated min")),

  pixels.sub %>%
  group_by(treatment) %>%
  summarise_at(vars(covTSF1pf:covTSF10, -Level3Ecoregions
                    ),max, na.rm = TRUE) %>%
  mutate(stat= ifelse(treatment==1, "treated max", "untreated max"))) %>%
  dplyr::select(-treatment) ))

colnames(summ.stats) <-c("untreated mean",
                        " treated mean",
                         "untreated sd",
                        "treated sd",
                        "untreated min",
                        "treated min",
                         "untreated max",
                       "treated max")

# Save summary stats
write.csv(summ.stats, "tables/summarystatisticsnew.csv")



### Format data for difference-in-differences model: #######
pixels.sub$fireID <- paste(pixels.sub$Fire_Year, pixels.sub$Shape_Area, sep="-")

time0 <- pixels.sub
time0$time <- 0
time0$sagecov <- time0$covTSF0
time0$pixelid <- as.factor(seq.int(nrow(time0)))

time10 <- pixels.sub
time10$time <- 1
time10$sagecov <- time10$covTSF10
time10$pixelid <- as.factor(seq.int(nrow(time10)))

unmatchdatDID <- rbind(time0, time10)
unmatchdatDID$DID <- unmatchdatDID$treatment*unmatchdatDID$time

# indices for hierarchical structure:
unmatchdatDID <-unmatchdatDID %>%
  transform(fireindex=as.numeric(factor(unmatchdatDID$fireID))) %>%
  transform(pixelindex=as.numeric(factor(unmatchdatDID$pixelid)))

pixels.in.fires <- unique(unmatchdatDID[,c('pixelindex','fireindex')])
colnames(pixels.in.fires)[2] <- "PixelinFire"

length(unique(pixels.sub$fireID))
length(unique(matchdat$fireID))


## Format data for panel regression:

pixels.sub2 <- pixels.sub
pixels.sub2$uniqueID <- seq.int(nrow(pixels.sub2))

sagecovers <- pixels.sub2 %>% 
  dplyr::select(treatment:covTSF10, 
                Treatment_status, uniqueID, fireID) %>%
  pivot_longer(cols=starts_with("covTSF"), 
               names_to = "yearTSF", values_to = "sagecov")%>%
  relocate(c(uniqueID, fireID, 
             Fire_Year, Treatment_status)) %>%
  mutate(yearTSF=replace(yearTSF, yearTSF=="covTSF1pf", 
                         "covTSF-1")) %>%
  mutate(yearTSF=str_sub(yearTSF, 7, 9)) %>%
  filter(yearTSF!=-1)

ppt <- pixels.sub2 %>% 
  dplyr::select(uniqueID, ppt0:ppt10) %>%
  pivot_longer(cols=starts_with("ppt"), 
               names_to = "yearTSF", values_to = "weather.sprppt") %>%
  mutate(yearTSF=str_sub(yearTSF, 4, 5))

temp <- pixels.sub2 %>% 
  dplyr::select(uniqueID, temp0:temp10) %>%
  pivot_longer(cols=starts_with("temp"), 
               names_to = "yearTSF", values_to = "weather.sprtmean") %>%
  mutate(yearTSF=str_sub(yearTSF, 5, 6))


paneldat <- sagecovers %>% left_join( ppt, 
                                      by=c("uniqueID", "yearTSF") )%>% 
  left_join( temp, by=c("uniqueID", "yearTSF") )

# Match ppt and lagged ppt values to appropriate years
ppt$yearTSF <- as.numeric(ppt$yearTSF)+1
ppt$lag.sprppt <- ppt$weather.sprppt
pptlag <- ppt %>% dplyr::select(uniqueID, yearTSF, lag.sprppt)

# Match temp values to appropriate years
temp$yearTSF <- as.numeric(temp$yearTSF)+1
temp$lag.sprtmean <- temp$weather.sprtmean
templag <- temp %>% dplyr::select(uniqueID, yearTSF, lag.sprtmean)

paneldat$yearTSF <- as.numeric(paneldat$yearTSF)

paneldat <- paneldat %>% left_join( pptlag, 
                                    by=c("uniqueID", "yearTSF") )%>% 
                          left_join( templag, 
                                     by=c("uniqueID", "yearTSF") )

paneldat$treatmentdummy <- ifelse(paneldat$yearTSF>0 & 
                                    paneldat$treatment==1,
                                  1, 0)

paneldat$fireID <- paste(paneldat$Fire_Year, 
                         paneldat$Shape_Area, sep="-")


# Save formatted datasets:
write.csv(pixels.sub, "formatteddata/randomsites_new.csv")
write.csv(unmatchdatDID, "formatteddata/unmatchDIDdata_new.csv")
write.csv(paneldat, "formatteddata/paneldata_new.csv")


rm(pixels.treat, pixels.untreat,data, data2, sage)

