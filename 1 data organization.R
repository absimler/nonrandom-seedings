### Nonrandom treatment locations and sagebrush seeding success #############
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

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
pixels.treat <- read.csv("formatteddata/pixels.treat.csv")
pixels.untreat <- read.csv("formatteddata/pixels.untreat.csv")


### b) Randomly select treated and untreated pixels: #######

### List of covariates that require values for analysis:
covars <- c('covTSF1pf',  "clay", "sand", 'febaprppt.pr',  
            'febaprtmean.tmean', 
            'novjanppt.pr', "novjantmean.tmean", 
            'Shape_Area', 'distroad', 'heatload',  
            "elevation", 
            'Level3Ecoregions', "Fire_Year",
            'covTSF0', "covTSF1", "covTSF2", "covTSF3", "covTSF4", "covTSF5",
            "covTSF6", "covTSF7", "covTSF8", "covTSF9", "covTSF10", "ppt1pf",
            'ppt0', "ppt1", "ppt2", "ppt3", "ppt4", "ppt5",
            "ppt6", "ppt7", "ppt8", "ppt9", "ppt10", "temp1pf",
            'temp0', "temp1", "temp2", "temp3", "temp4", "temp5",
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


### c) Data summary statistics #######

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



### d) Format data for difference-in-differences model: #######
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


##e) Format data for panel regression: #####

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

ppt$yearTSF <- as.numeric(ppt$yearTSF)+1
ppt$lag.sprppt <- ppt$weather.sprppt
pptlag <- ppt %>% dplyr::select(uniqueID, yearTSF, lag.sprppt)

temp$yearTSF <- as.numeric(temp$yearTSF)+1
temp$lag.sprtmean <- temp$weather.sprtmean
templag <- temp %>% dplyr::select(uniqueID, yearTSF, lag.sprtmean)


paneldat$yearTSF <- as.numeric(paneldat$yearTSF)

paneldat <- paneldat %>% left_join( pptlag, 
                                    by=c("uniqueID", "yearTSF") )%>% 
  left_join( templag, by=c("uniqueID", "yearTSF") )

paneldat$treatmentdummy <- ifelse(paneldat$yearTSF>0 & 
                                    paneldat$treatment==1,
                                  1, 0)

paneldat$fireID <- paste(paneldat$Fire_Year, paneldat$Shape_Area, sep="-")



# f) Save formatted datasets for model fits ######
write.csv(pixels.sub, "formatteddata/randomsites_final.csv")
write.csv(unmatchdatDID, "formatteddata/unmatchDIDdata_final.csv")
write.csv(paneldat, "formatteddata/paneldata_final.csv")


rm(pixels.treat, pixels.untreat, ppt, temp,templag, time0, time10, summ.stats,
   pptlag,pixels.sub2)

