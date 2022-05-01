### Nonrandom treatment locations and sagebrush seeding success##########
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

#### Supplementary: Cluster robust standard errors comparison #####

library(MASS)
library(sandwich)
library(coefplot)
library(lmtest)
library(broom)
library(fixest)

#### cluster robust standard errors, using glm(), fenbin() and sandwich() ####

### DID with CRSE####
clustrobust.did <- glm.nb(sagecov~treatment +
                                    time + DID, 
                                  data=diddat)
summary(clustrobust.did)
coefplot(clustrobust.did)

coeftest(clustrobust.did, 
         vcov = vcovCL, 
         cluster = ~ pixelid + fireID)

clusterest <- tidy(coeftest(clustrobust.did, 
                       vcov = vcovCL, 
                       cluster = ~ pixelid + fireID),
                   conf.int=T, conf.level=0.95)

write.csv(clusterest, "tables/clusterrobustDiD.csv")




## WE Panel model w/CRSE ###
clustrobust.wep <- glm.nb(sagecov~treatment + as.factor(yearTSF) + 
                            treatmentdummy + 
                            scale(weather.sprppt) +
                            scale(weather.sprtmean) + as.factor(fireID), 
                          data=paneldat )


coeftest(clustrobust.wep, 
         vcov = vcovCL, 
         cluster = ~ uniqueID + fireID)

clusterest.pan  <- tidy(coeftest(clustrobust.wep, 
                            vcov = vcovCL, 
                            cluster = ~ uniqueID + fireID),
                   conf.int=T, conf.level=0.95)


write.csv(clusterest.pan, "tables/clusterrobustpanel.csv")




#### Compare predicted average treatment effects for treated group from CRSE v. varying intercept models #####

didcrse.yhattreat <- predict(clustrobust.did, 
                             newdata=data.frame(DID = 1, treatment=1,
                              time = 1) , se.fit=T, type="response")
didcrse.yhatuntreat <- predict(clustrobust.did,
                               newdata=data.frame(DID = 0, treatment=1,
                                                  time = 1), se.fit=T, type="response")


didcrse.yhattreat$fit - didcrse.yhatuntreat$fit


didcrse.yhattreat
didcrse.yhatuntreat




wepcrse.yhattreat <- predict(clustrobust.wep, 
                             newdata=data.frame(treatmentdummy = 1, treatment=1,
                                                yearTSF = 10,
                                                weather.sprppt=mean(paneldat$weather.sprppt),
                                                weather.sprtmean=mean(paneldat$weather.sprtmean)), 
                             se.fit=T, type="response")

wepcrse.yhatuntreat <- predict(clustrobust.wep,
                               newdata=data.frame(treatmentdummy = 0, treatment=1,
                                                  yearTSF = 10,
                                                  weather.sprppt=mean(paneldat$weather.sprppt),
                                                  weather.sprtmean=mean(paneldat$weather.sprtmean)),
                               se.fit=T, type="response")


wepcrse.yhattreat$fit - wepcrse.yhatuntreat$fit



