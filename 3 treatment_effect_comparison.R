### Nonrandom treatment locations and sagebrush seeding success ########
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######


#### PART 3) Compare restoration treatment effects identified by statistical approaches ###

library(rstan)
library(rstanarm)
library(brms)
library(tidyverse)
library(cmdstanr)
library(performance)

options(mc.cores = parallel::detectCores())

## Load formatted datasets ####
# Same initial set of randomly selected sagebrush sites, before or after matching, with or without panel structure

sagesites <- read.csv("formatteddata/randomsites_new.csv")
matchdat <- read.csv("formatteddata/matcheddata_new.csv")
diddat <- read.csv("formatteddata/unmatchDIDdata_new.csv")
paneldat <- read.csv("formatteddata/paneldata_new.csv")

length(unique(diddat$fireID)) # sample sizes for indiv. fires
length(unique(matchdat$fireID))

###### COMPARISON OF MODELS' ESTIMATED TREATMENT EFFECTS #######
### "Naive" null model  ####
## Treatment effect estimated from simple difference in means, without matching or covariates 

treatmenteff.nomatch <- brm(covTSF10~treatment, 
                                    data=sagesites, 
                            chains=4, iter=2000, 
                            family="negbinomial"(link="log"))

saveRDS(treatmenteff.nomatch, 
        file = "modelfits/treatmenteff.nomatch_final.rda")


### "Covariate control model" ##
# Treatment effect estimated with covariates hypothesized to influence treatment groups

treatmenteff.covars <- brm(covTSF10~treatment+
                                 scale(wintersprppt) + 
                                 scale(I(wintersprppt)^2) + 
                                 scale(febaprtmean.tmean) + 
                                 scale(I(febaprtmean.tmean)^2)+
                                 scale(clay) + scale(sand)+
                                 scale(heatload) +
                                 scale(I(heatload^2)) +
                                 scale(elevation) +
                                 scale(I(elevation^2)) +
                                 Level3Ecoregions + (1 | fireID), 
                                 data=sagesites, 
                                 chains=4, iter=2000,
                                 family="negbinomial"(link="log"),
                           backend = "cmdstanr", 
                           threads = threading(8))

saveRDS(treatmenteff.covars, 
        file = "modelfits/treatmenteff.covarsfireid_final.rda")

mcmc_plot(treatmenteff.covars, 
          regex = "^b", prob_outer=0.95)


#### "Matched model" ###
# Treatment effect estimated following PSM: ####
treatmenteff.match <- brm(covTSF10~treatment, 
                          data=matchdat, 
                          chains=4, iter=2000,
                          family="negbinomial"(link="log"),
                          backend = "cmdstanr", 
                          threads = threading(8))

saveRDS(treatmenteff.match, 
        file = "modelfits/treatmenteff.match_final.rda")


## Difference in differences model ######
## Incorporates hierarchical structure for site and fire identity to account for clustering of observations:

treatmenteff.did <- brm(sagecov~treatment + time +
                              DID + (1|pixelid) +
                              (1|fireID), 
                              data=diddat,
                              chains=4,
                              iter=2000,
                              family="negbinomial"(link="log"),
                              backend = "cmdstanr", 
                              threads = threading(8))

saveRDS(treatmenteff.did, 
        file = "modelfits/treatmenteff.did.sitefireid_final.rda")


### Within-estimator panel regression ####
paneldat$treatmentdummy <- as.factor(paneldat$treatmentdummy)
paneldat$treatment <- as.factor(paneldat$treatment)
paneldat$yearTSF <- as.factor(paneldat$yearTSF)

## model with Fire Identity
treatmenteff.panel <- brm(sagecov~
                                 treatment+
                     treatmentdummy+
                     yearTSF +
                     scale(weather.sprppt) + 
                     scale(weather.sprtmean) +
                     (1|uniqueID) + 
                      (1|fireID),
                   data=paneldat, chains=4, init=0, thin=2,
                   iter=2000, family="negbinomial"(link="log"),
                   backend = "cmdstanr", 
                   threads = threading(8))

saveRDS(treatmenteff.panel,
        file = "modelfits/treatmenteff.panelsitefireid_final.rda")


### Efficacy of treatments in different environmental contexts ####

treateff.did.environ <- brm(sagecov ~ 
                         treatment + 
                         time + 
                         DID +
                         DID*scale(wintersprppt) + 
                         DID*scale(I(wintersprppt^2)) +
                         DID*scale(febaprtmean.tmean) +
                         DID*scale(I(febaprtmean.tmean^2))+
                         DID*scale(clay) +
                         DID*scale(I(clay^2)) + 
                           DID*scale(sand) +
                           DID*scale(I(sand^2)) +
                           (1|Level3Ecoregions) +
                           (1|fireID)+
                        (1|pixelid),
                       data=diddat, 
                      chains=4, 
                      iter=2000, thin=2,
                      family="negbinomial"(link="log"),
                      backend = "cmdstanr", 
                   threads = threading(8))

saveRDS(treateff.did.environ, 
    file = "modelfits/treatmenteff.did.environ_fireidpixelid.rda")


### Fit metrics:
mae(treatmenteff.nomatch)
mae(treatmenteff.match)
mae(treatmenteff.covars)
mae(treatmenteff.did)
mae(treatmenteff.panel.small)


pp_check(treatmenteff.panel)
pp_check(treatmenteff.did)

