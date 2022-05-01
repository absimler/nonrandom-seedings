### Nonrandom treatment locations and sagebrush seeding success #############
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######


#### Figures for comparing treatment effects among statistical approaches ####

library(rstan)
library(rstanarm)
library(tidybayes)
library(brms)
library(tidyverse)
library(bayesplot)
library(cowplot)
library(modelr)

# Load model fits & data:
treatmenteff.nomatch <- readRDS("modelfits/treatmenteff.nomatch_feb28.rda")
treatmenteff.covars <- readRDS("modelfits/treatmenteff.covarsfireid_feb28.rda")
treatmenteff.match <- readRDS("modelfits/treatmenteff.match_feb28.rda")
treatmenteff.did <- readRDS( "modelfits/treatmenteff.did.sitefireid_feb28.rda")
treatmenteff.panel <- readRDS( "modelfits/treatmenteff.panelsitefireid_feb28.rda")
did.env <- readRDS("modelfits/treatmenteff.did.environ_mar2.rda")

sagesites <- read.csv("formatteddata/randomsites_new.csv")
DIDdat <- read.csv("formatteddata/unmatchDIDdata_new.csv")
matchdat <- read.csv("formatteddata/matcheddata_new.csv")
paneldat <- read.csv("formatteddata/paneldata_new.csv")


##### FIGURE: Comparison of treatment effects at treated sites ################

## Treatment effect for nonlinear model:
## = exp(time + group + interaction + intercept/covariates) - exp(time + treatment + covariates/intercept)


#Unmatched model - posterior samples turned into treat effs:
modfit.nomatch <- posterior_samples(treatmenteff.nomatch, "^b")

unmatch.treat <- exp(modfit.nomatch[1] + modfit.nomatch[2])
unmatch.treat$untreated <- exp(modfit.nomatch[1])

unmatch.treat$treateff<- unmatch.treat$b_Intercept - unmatch.treat$untreated
unmatch.treat$model <- "Naive model\n(Difference in means)"
unmatch.treat <- unmatch.treat[c("treateff","model")]

#Env covariate model
modfit.env <- data.frame(treatmenteff.covars) 
modfit.env <-modfit.env %>% dplyr::select(b_Intercept, b_treatment)

covars.treat <- exp(modfit.env[1] + modfit.env[2])
covars.treat$untreated <- exp(modfit.env[1] )

covars.treat$treateff<- covars.treat$b_Intercept - 
  covars.treat$untreated
covars.treat$model <- "Regression with\nenvironmental covariates\n(with fire ID)"
covars.treat <- covars.treat[c("treateff","model")]


#Matched model
modfit.match <- posterior_samples(treatmenteff.match) 

match.treat <- exp(modfit.match[1] + modfit.match[2])
match.treat$untreated <- exp(modfit.match[1])

match.treat$treateff<- match.treat$b_Intercept - match.treat$untreated
match.treat$model <- "Regression following\npropensity score matching"
match.treat <- match.treat[c("treateff","model")]


## Treatment effect for DID model
modfit.did <- posterior_samples(treatmenteff.did, "^b") 

did.treat <- exp(modfit.did[1] + 
                   modfit.did[2] + 
                   modfit.did[3] +
                   modfit.did[4])
did.treat$untreated <- exp(modfit.did[1] +
                             modfit.did[2] + 
                             modfit.did[3])

did.treat$treateff<-did.treat$b_Intercept-did.treat$untreated
did.treat$model <- "Difference-in-Differences\n(with location & fire ID)"
did.treat <- did.treat[c("treateff","model")]


## Treatment effect for panel model:
modfit.panel <- posterior_samples(treatmenteff.panel, "^b") 

panel.treat <- exp(modfit.panel[1] + 
                     modfit.panel[2] + 
                     modfit.panel[3] +
                   modfit.panel[13])
panel.treat$untreated <- exp(modfit.panel[1] +
                             modfit.panel[2] + 
                               modfit.panel[13])

panel.treat$treateff<- panel.treat$b_Intercept - 
  panel.treat$untreated
panel.treat$model <- "Within-estimator\npanel regression\n(with location & fire ID)"
panel.treat <- panel.treat[c("treateff","model")]


# Combine all predicted ATEs across models:
ate.preds <- rbind(panel.treat,
                   did.treat, 
                   unmatch.treat, 
                   covars.treat,
                   match.treat)

ate.preds <- data.frame(ate.preds$treateff$b_Intercept, 
                        ate.preds$model)
colnames(ate.preds) <- c("treateff", "model")

ate.preds$model <- factor(ate.preds$model, 
          levels = c("Within-estimator\npanel regression\n(with location & fire ID)", 
                                                   "Difference-in-Differences\n(with location & fire ID)", 
                                                   "Regression with\nenvironmental covariates\n(with fire ID)",
                                                   "Regression following\npropensity score matching",
                                                   "Naive model\n(Difference in means)"))


treateff.pred <-ggplot(ate.preds, 
                aes(x=as.factor(model),
                        y=treateff, color=model, 
                        fill=model))+
                geom_hline(yintercept = 0, 
                           linetype = "dashed", 
                           size=1, color = "grey50") +
               stat_halfeye(fill="grey80", 
                            .width = c(0.5, 0.95)) +
                scale_fill_manual(name="")+
               scale_color_manual(name="",
                           values = c(  "#30156d","#30156d",
                         "grey30", "grey30", "#CC3333")) +
                      xlab("Model Type") + 
                      ylab("Average treatment effect for treated sites\n(Change in sagebrush percent cover, 10 years post-treatment)") +
                   theme(legend.position="none", 
                         legend.box = "horizontal") + 
                      theme_bw() +
  theme(legend.title = element_blank(),
        legend.position=c(0.15, 0.85),
        axis.text.y = element_text( size=8), 
        axis.text.x= element_text(size=8),
        axis.title = element_text(size=8),
        legend.text = element_text(size=8)) + 
    guides(fill=FALSE, color=FALSE) + 
  coord_flip()

treateff.pred


comboATT <- plot_grid(treateff.pred, paneltime, 
                      ncol=1, nrow=2, labels=c("a", "b"))

comboATT


att.summ <-  as.data.frame(bind_rows(
                 ate.preds %>%
                    group_by(model) %>%
                    summarise_at(vars(treateff),
                                 median, na.rm = TRUE) %>%
                   mutate(stat= "median"),
                 ate.preds %>%
                   group_by(model) %>%
                   summarise_at(vars(treateff),
                                mean, na.rm = TRUE) %>%
                   mutate(stat= "mean"),
                     ate.preds %>%
                       group_by(model) %>%
                       summarise_at(vars(treateff),
                                    sd, na.rm = TRUE) %>%
                   mutate(stat= "sd")))

write.csv(att.summ, "tables/avgtreatoftreated.csv")



##### Parameter plots Environmental model #################

posteriorall <- mcmc_intervals_data(did.environ, 
                                    prob_outer=0.95,
                                 prob=0.5)
posterior <- posteriorall%>% 
              filter(grepl('b_', parameter)) %>%
                mutate(reorder = c(29:1))
                       
                       
posterior$nonzero <- NA
posterior$nonzero[posterior$ll>0 & posterior$hh>0] <- "nonzero"
posterior$nonzero[posterior$ll<0 & posterior$hh<0] <- "nonzero"
posterior$nonzero[is.na(posterior$nonzero)] <- "zero"



varnames1 <- c("Intercept", "Group indicator", 
               "Time indicator", "Treatment effect", 
               "Nov-Apr precipitation",
               "Nov-Apr precipitation\n (quadratic)", 
               "Feb-Apr mean\ntemperature",
               "Feb-Apr mean\ntemperature (quadratic)",
               "soil percent clay",
               "soil percent\nclay (quadratic)",
               "Central Basin\nand Range ER", 
               "Columbia Plateau ER",
               "Eastern Cascades\nFoothills ER",
               "High Plains ER",
               "Idaho Batholith ER",
               "Middle Rockies ER",
               "Northern Basin\nand Range ER", 
               "NW Great Plains ER",
               "Sierra Nevada ER",
               "Snake River\nPlain ER", 
               "Southern Rockies ER",
               "Wasatch/Uinta\nMountains ER",
               "Wyoming Basin ER", 
               "Treatment*Nov-Apr\nprecipitation", 
               "Treatment*Nov-Apr\nprecipitation (quadratic)",
               "Treatment*Feb-Apr\ntemperature", 
               "Treatment*Feb-Apr\ntemperature (quadratic)",
               "Treatment*soil percent clay",
               "Treatment*soil percent\nclay (quadratic)")
posterior$varnames <- c(varnames1)

#posterior1 <- filter(posterior, grepl('Treatment*', varnames))
posterior1 <- filter(posterior, !grepl('Ecoregions', parameter))
posterior2 <- filter(posterior, grepl('Ecoregions', parameter))

# posterior parameter estimate plots
postplot.env <-  ggplot(subset(posterior), 
                    aes(x = reorder(varnames, reorder),
                        linetype=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  shape = 16, size = 3/4) +
  scale_linetype_manual(values=c("solid", "dashed"))+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y = element_text( size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  theme(legend.title=element_blank(), 
        legend.position=c(0.75, 0.75)) +
  xlab(NULL) +
  ylab("Effect on sagebrush cover")+
  guides(linetype=FALSE) 

postplot.env



postplot2.env <-  ggplot(posterior1, 
                         aes(x = reorder(varnames, 
                                         reorder),
                                     linetype=nonzero)) +
  geom_hline(yintercept = 0, 
             linetype = 3, size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  shape = 16, size = 3/4) +
  scale_linetype_manual(values=c("solid", "dashed"))+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y = element_text( size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  theme(legend.position=c(0.85, 0.75), 
        legend.title=element_blank()) +
  xlab(NULL) +
  ylab("Effect of parameter on probability\nof receiving seeding")+
  guides(linetype=FALSE, color=FALSE) 

postplot2.env


##### 10) Environment*Treatment Marginal effects plots #######
mytheme <-  theme_bw() + 
  theme(axis.text.y = element_text( size=8), 
        axis.text.x= element_text(size=8),
        axis.title = element_text(size=8), 
        legend.text=element_text(size=8))


## precip effects:
ppt.gradient <- rep(seq(from=min(sagesites$wintersprppt, 
                                 na.rm=T), 
                        to= max(sagesites$wintersprppt, 
                                na.rm=T),
                        length.out=500), 2)

preds.ppt <- DIDdat %>%
  data_grid(time = 1,
            treatment = 1,
            DID = c(rep(1, 100), 
                    rep(0, 100)),
            wintersprppt = ppt.gradient,
            febaprtmean.tmean = mean(febaprtmean.tmean, 
                                     na.rm=T),
            clay = mean(clay, na.rm=T),
            Level3Ecoregions="Central Basin\nand Range") %>%
  add_epred_draws(did.environ, draws=100, re_formula=NA)


cf.ppt.recovery <-   ggplot(preds.ppt, 
                            aes(x = wintersprppt,
                            y = .epred, 
                            fill=as.factor(DID), 
                            color=as.factor(DID))) +
  stat_lineribbon(.width = c(.5), alpha = 0.5, size=1.5) +
  stat_lineribbon(.width = c(.0001), size=1.5) +
  scale_color_manual(values=c( "#0066CC", "#CC3333"), 
            labels=c("E[if untreated]", "E[if treated]"))  +
  scale_fill_manual(values=c( "#0066CC", "#CC3333"),
            labels=c("E[if untreated]", "E[if treated]"))  +
  ylab("Sagebrush cover \n 10 years after seeding")+ 
  xlab("November-April total precipitation (mm)") +
  mytheme + theme(legend.title=element_blank(), 
                  legend.position = c(0.2, 0.85)) +
  theme(rect = element_rect(fill = "transparent"))

cf.ppt.recovery


preds.treat <- subset(preds.ppt, DID==1)
preds.untreat <- subset(preds.ppt, DID==0)

preds.treat$.value.untreat <- preds.untreat$.epred
preds.treat$cover.diff <- preds.treat$.epred - preds.treat$.value.untreat

probtreatbefore <- readRDS(file = "modelfits/treatmentprob_beforePSM_feb28.rda")

treatpreds <- sagesites %>%
  data_grid(wintersprppt = ppt.gradient,
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = mean(febaprtmean.tmean, 
                                     na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions = "Central Basin\nand Range") %>%
  add_fitted_draws(probtreatbefore)

avgprediction <-  treatpreds %>% group_by(wintersprppt) %>%
  median_qi(.value) %>% filter(.value<0.51 & .value>0.49)

cf.ppt.compare <-   ggplot(preds.treat, 
                    aes(x = wintersprppt, y = cover.diff,
                        fill = stat(abs(x) <= 911),
                        color = stat(abs(x) <= 911 ))) +
  stat_lineribbon(.width = c(.5), alpha = 0.5, size=2) +
  stat_lineribbon(.width = c(.00001), size=2) +
    scale_fill_manual(values=c( "#0066CC", "#CC3333"), 
          labels=c("Location unlikely to be treated (<50%)",     "Location likely to be treated (>50%)"))  +
    scale_color_manual(values=c( "#0066CC", "#CC3333"), 
       labels=c("Location unlikely to be treated (<50%)", 
          "Location likely to be treated (>50%)"))  +
  ylab("Predicted difference in sagebrush cover\nassociated with seeding")+ 
  xlab("November-April total precipitation (mm)") +
  geom_hline(yintercept=0)+
  mytheme + theme(legend.title=element_blank(), 
                legend.position = c(0.49, 0.18)) +
  theme(rect = element_rect(fill = "transparent"))
cf.ppt.compare

avgprediction <-  preds.treat %>% 
  group_by(wintersprppt) %>%
  median_qi(cover.diff) 

avgprediction[which.max(avgprediction$cover.diff),]
avgprediction[which.min(avgprediction$cover.diff),]


## temp effects:
temp.gradient <- rep(
  seq(from=min(sagesites$febaprtmean.tmean,na.rm=T), 
      to= max(sagesites$febaprtmean.tmean, na.rm=T),
                         length.out=500), 2)


preds.temp <- DIDdat %>%
  data_grid(time = 1,
            treatment = 1,
            DID = c(rep(1, 100), rep(0, 100)),
            wintersprppt = mean(wintersprppt, na.rm=T),
            febaprtmean.tmean = temp.gradient,
            distroad = mean(distroad, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            clay = mean(clay, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            Level3Ecoregions = "Central Basin\nand Range") %>%
  add_fitted_draws(did.environ, draws=100, re_formula=NA)


cf.temp.recovery <-   ggplot(preds.temp, aes(x = febaprtmean.tmean, y = .value, fill=as.factor(DID), color=as.factor(DID))) +
  stat_lineribbon(.width = c(.5), alpha = 0.5, size=1.5) + 
  stat_lineribbon(.width = c(.0001), size=1.5) + 
  scale_color_manual(values=c( "#197bdd", "#CC3333"), 
                     labels=c("E[if untreated]", "E[if treated]"))  +
  scale_fill_manual(values=c( "#197bdd", "#CC3333"),
                    labels=c("E[if untreated]", "E[if treated]"))  +
  ylab("Sagebrush cover \n 10 years after seeding")+ 
  xlab("February-April mean temperature (C)") +
  #xlim(-5,8) +
  mytheme + theme(legend.title=element_blank(), legend.position = c(0.7, 0.7)) +
  theme(rect = element_rect(fill = "transparent"))

cf.temp.recovery


preds.treat <- subset(preds.temp, DID==1)
preds.untreat <- subset(preds.temp, DID==0)

preds.treat$.value.untreat <- preds.untreat$.value
preds.treat$cover.diff <- preds.treat$.value - preds.treat$.value.untreat

treatpreds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = temp.gradient,
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions = "Central Basin\nand Range") %>%
  add_fitted_draws(probtreatbefore)


avgprediction <-  treatpreds %>% group_by(febaprtmean.tmean) %>%
  median_qi(.value) %>% filter(.value<0.505 & .value>0.495)


cf.temp.compare <-   ggplot(preds.treat, 
                            aes(x = febaprtmean.tmean, 
                             y = cover.diff,
                             fill = stat(x > 0.59),
                             color = stat(x> 0.59))) +
  stat_lineribbon(.width = c(.5), alpha = 0.5, size=2) +
  stat_lineribbon(.width = c(.0001), size=2) +
  scale_fill_manual(values=c( "#0066CC", "#CC3333"), 
         labels=c("Location unlikely to be treated (<50%)", 
         "Location likely to be treated (>50%)"))  +
  scale_color_manual(values=c( "#0066CC", "#CC3333"), 
         labels=c("Location unlikely to be treated (<50%)", 
            "Location likely to be treated (>50%)"))  +
  ylab("Predicted difference in sagebrush cover\nassociated with seeding")+ xlab("February-April mean temperature (C)") +
  geom_hline(yintercept=0)+
  mytheme + theme(legend.title=element_blank(), 
                  legend.position = c(0.45, 0.2))+
  theme(rect = element_rect(fill = "transparent"))
cf.temp.compare

avgprediction <-  preds.treat %>% group_by(febaprtmean.tmean) %>%
  median_qi(cover.diff) 

avgprediction[which.max(avgprediction$cover.diff),]
avgprediction[which.min(avgprediction$cover.diff),]


### 11) Combined plots#####

precips <- plot_grid(cf.ppt.recovery, 
                     cf.ppt.compare, ncol=2, nrow=1, 
                     labels=c("a", "b"))
temps <- plot_grid(cf.temp.recovery, cf.temp.compare, 
                   ncol=2, nrow=1, labels=c("c", "d"))

cfs <- plot_grid(precips, temps, ncol=1, nrow=2)

envplots <- plot_grid(postplot.env, cfs, 
            ncol=1, nrow=2, rel_heights =c(1.5, 2),
                      labels=c("A", ""))
envplots



### Save all plots #####

pdf('figures/environmental var final.pdf',
     width=7, height=7
)
cfs
dev.off()


jpeg('figures/environmental var pars final.jpeg',
     width=9, height=11, units="in", res=600
)
postplot.env
dev.off()


pdf('figures/combo att final.pdf',
    width=6, height=7
)
comboATT
dev.off()

