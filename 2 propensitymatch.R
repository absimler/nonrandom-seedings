### Nonrandom treatment locations and sagebrush seeding success #############
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

#########################################################

#### PART 2) Propensity scores and comparison of biophysical characteristics between treated and untreated pixels ######


library(MatchIt)
library(rstan)
library(rstanarm)
library(brms)
library(tidyverse)
library(cmdstanr)

sagesites <- read.csv("formatteddata/randomsites.csv")
sagesites %>% count(Level3Ecoregions)

### a) Assess potential drivers of probability of a pixel of receiving treatment: ####

sagesites$Level3Ecoregions <- as.factor(sagesites$Level3Ecoregions)

probtreatbefore <- brm(treatment~ 
                         scale(covTSF1pf) + scale(I(covTSF1pf^2)) +
                         scale(covTSF0) +
                         scale(wintersprppt) + 
                         scale(I(wintersprppt^2)) +
                         scale(febaprtmean.tmean) + 
                         scale(I(febaprtmean.tmean^2)) + 
                         scale(Shape_Area) +
                         scale(clay) + scale(sand)+
                         scale(elevation) +
                         scale(I(elevation^2)) +
                         scale(distroad) +
                         scale(I(distroad^2)) + 
                         (1|Level3Ecoregions),
                       family=bernoulli(link=logit), 
                       data=sagesites,
                       chains=2, iter=2000,  
                       control=list(max_treedepth=30),
                       backend = "cmdstanr", 
                       threads = threading(8))

plot(probtreatbefore, prob_outer=0.95) 

## Bayesian R-squared for binomial glm
proboftreat.R2 <- median(bayes_R2(probtreatbefore))
proboftreat.R2

saveRDS(probtreatbefore, 
        file = "modelfits/treatmentprob_beforePSM.rda")
  


### b) Complete Propensity Score Matching: ######

#Set caliper to 0.25*Std deviation of the propensity scores:
probtreatbefore <- readRDS(file = 
                  "modelfits/treatmentprob_beforePSM.rda")

caliper = 0.2*sd(colMeans(posterior_linpred(probtreatbefore, 
                                             transform = TRUE)))

psmmodel <- matchit(treatment ~ covTSF1pf + 
                        I(covTSF1pf^2) +
                        covTSF0 +
                        wintersprppt + 
                        I(wintersprppt^2) +
                        febaprtmean.tmean + 
                        I(febaprtmean.tmean^2) + 
                        novjanppt.pr + 
                        I(novjanppt.pr^2) +
                        novjantmean.tmean + 
                        I(novjantmean.tmean^2) + 
                        elevation + 
                        I(elevation^2)+ 
                        heatload + 
                        I(heatload^2) +
                        Shape_Area + 
                        I(distroad^2) + 
                        clay + sand +
                        distroad + 
                        Level3Ecoregions,
                      method = "nearest", discard="both", 
                    caliper=caliper, data = sagesites)

matches <- match.data(psmmodel) ## Subset data to just those that matched


### d) Assess effects after matching using glm: ####
probtreatafter <- brm(treatment~ 
                        scale(covTSF1pf) + 
                        scale(I(covTSF1pf^2)) +
                        scale(covTSF0) +
                        scale(wintersprppt) + 
                        scale(I(wintersprppt^2)) +
                        scale(febaprtmean.tmean) + 
                        scale(I(febaprtmean.tmean^2)) +
                        scale(Shape_Area) + 
                        scale(clay) + 
                        scale(sand)+
                        scale(elevation) +
                        scale(I(elevation^2)) +
                        scale(distroad) +
                        scale(I(distroad^2)) +
                        (1|Level3Ecoregions),
                      family=bernoulli(link=logit), 
                      data=matches,
                      chains=4, iter=2000, cores=4,  
                      control=list(max_treedepth=30, adapt_delta=0.99),
                      backend = "cmdstanr", 
                      threads = threading(8))

plot(probtreatafter, prob_outer=0.95)

saveRDS(probtreatafter, 
        file = "modelfits/treatmentprob_afterPSM.rda")

probtreatafter <- readRDS(file = 
                             "modelfits/treatmentprob_afterPSM_feb28.rda")

### e) Match sets and save matched data #######

matches$included <- TRUE
matchdat <- merge(sagesites, matches, all=T)
matchdat <- subset(matchdat, included==TRUE)
matchdat$fireID <- paste(matchdat$Fire_Year, matchdat$Shape_Area)

write.csv(matchdat, "formatteddata/matcheddata_new.csv")
