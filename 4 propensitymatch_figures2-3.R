### Nonrandom treatment locations and sagebrush seeding success #############
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

## PART 3) Figures associated with Propensity Score analysis and drivers of selection bias:#######

matchdat <- read.csv("formatteddata/matcheddata_final.csv")

mytheme <-theme(axis.text.y = 
                  element_text( size=7),
        axis.title = element_text(size=7), 
        axis.text.x=element_text(size=7))


## a) Compare propensity scores for matched and unmatched pixels: ####
propensity.unmatched <- data.frame(pr_score = 
                           colMeans(posterior_linpred(probtreatbefore, 
                                                      transform = TRUE)),
                                   treatment = sagesites$treatment)

propensity.unmatched$Treatment_status <- ifelse(propensity.unmatched$treatment==1,
                                                "Treated", "Untreated")


unmatched.prop <- ggplot(propensity.unmatched, 
                         aes(x=pr_score, color=Treatment_status)) +
                          geom_density(aes(fill=Treatment_status), alpha=0.4) + 
                          xlab("Probability of receiving seeding") + 
                          scale_color_manual(values=c("#CC3333", "#197bdd")) +
                          scale_fill_manual(values=c("#CC3333", "#197bdd")) +
                          ggtitle("Before propensity matching") +
                          theme_bw() +   mytheme + 
                          theme(legend.title=element_blank(),  
                                             legend.position=c(0.25, 0.75),
                                             legend.text = element_text(size=7),
                                             plot.title = 
                                  element_text(hjust = 0.5, size=8)) + 
                             guides(fill=FALSE, color=FALSE) 

unmatched.prop

legend.prop <- get_legend(unmatched.prop)


propensity.matched <- data.frame(pr_score = 
          colMeans(posterior_linpred(probtreatafter, 
                                     transform = TRUE)),
                  treatment = matchdat$treatment)


matched.prop <- ggplot(propensity.matched, 
                       aes(x=pr_score, color=as.factor(treatment))) +
                        geom_density(aes(fill=as.factor(treatment)), alpha=0.4) +
                        scale_color_manual(values=c("#CC3333", "#197bdd"),
                                           labels=c("Untreated", "Treated")) +
                        scale_fill_manual(values=c("#CC3333", "#197bdd"),
                                          labels=c("Untreated", "Treated")) +
                        xlab("Probability of receiving seeding")+ 
                        ggtitle("After propensity matching") +
                        theme_bw() +mytheme+ 
                            theme(legend.title=element_blank(),    
                                  legend.text = element_text(size=7),
                                           legend.position=c(0.79, 0.75),
                                plot.title = element_text(hjust = 0.5, 
                                                          size=7))  


matched.prop


write.csv(propensity.matched, "sourcefigs/fig2a1.csv")
write.csv(propensity.unmatched, "sourcefigs/fig2a2.csv")

#Combined figure of density plots for propensity scores before and after PSM:
prop.grid <- plot_grid(unmatched.prop, matched.prop,
                       nrow=1, ncol=2, rel_widths = c(1, 1))
prop.grid



##b) Parameter estimates from propensity models before and after matching:#######

posterior <- bind_rows(mcmc_intervals_data(probtreatbefore, 
                                           prob_outer=0.95,
                                         prob=0.5) %>%
                                         mutate(reorder = c(46:1), 
                                         model="Unmatched\ndataset"),

                       mcmc_intervals_data(probtreatafter, 
                                           prob_outer=0.95,
                                       prob=0.5) %>%
                                       mutate(reorder = c(28:1), 
                                       model="Matched\nsubset") )

posterior$nonzero <- NA
posterior$nonzero[posterior$ll>0 & posterior$hh>0] <- "nonzero"
posterior$nonzero[posterior$ll<0 & posterior$hh<0] <- "nonzero"
posterior$nonzero[is.na(posterior$nonzero)] <- "zero"


posterior1 <- filter(posterior, 
                     grepl('scale', parameter) )
posterior2 <- filter(posterior, grepl('r_', parameter) | 
                       grepl('b_Intercept', parameter) )

varnames1 <- c( "Pre-fire sagebrush\ncover - linear", 
               "Pre-fire sagebrush\ncover - quadratic",
               "Post-fire sagebrush\ncover", 
               "Feb-Apr precipitation\n - linear",
               "Feb-Apr precipitation\n - quadratic", 
               "Feb-Apr mean\ntemperature - linear",
               "Feb-Apr mean\ntemperature - quadratic", 
               "Fire size", 
               "Soil percent clay", 
               "Soil percent sand",
               "Elevation", 
               "Elevation - quadratic",
               "Distance from a\nmajor road - linear",
               "Distance from a\nmajor road - quadratic")
posterior1$varnames <- c(varnames1, varnames1)


posterior2 <- posterior2 %>%
               mutate(varnames=str_sub(parameter,19, 35)) 

posterior2$varnames[posterior2$varnames==""] <- "Intercept"

# posterior parameter estimate plots
postplot <-  ggplot(posterior1, aes(x = reorder(varnames, reorder),
                                    color=model,
                                  shape=nonzero)) +
                geom_hline(yintercept = 0, linetype = 3, 
                           size=1, color = "#b0b5b3") +
                geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                                position= position_dodge(width=0.75),
                                size = 3/4) +
                scale_color_manual(name="",
                                   values = c("grey60", "#484c8d")) +
                scale_shape_manual(values=c(16, 17), 
                                   labels=c("95% CI does\nnot contain zero", 
                                            "95% CI\ncontains zero"))+
                
                coord_flip() +
                theme_bw() + 
                theme(axis.text.y = element_text( size=7), 
                      axis.text.x=element_text(size=7),
                      axis.title = element_text(size=7), 
                      legend.text = element_text(size=7)) +
                theme(legend.title=element_blank(), 
                      legend.position=c(0.79, 0.35)) +
                xlab(NULL) +
                ylab("Effect on probability\nof receiving seeding")+
                guides(linetype=FALSE) 

postplot



postplot2 <-  ggplot(posterior2, aes(x = reorder(varnames, reorder),
                                     color=model,
                                     shape=nonzero)) +
              geom_hline(yintercept = 0, linetype = 3, 
                         size=1, color = "#b0b5b3") +
              geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                              position= position_dodge(width=0.75),
                               size = 3/4) +
              scale_color_manual(name="",
                                 values = c("grey60", "#484c8d")) +
              scale_shape_manual(values=c(16, 17), 
                                 labels=c("95% CI does\nnot contain zero", 
                                          "95% CI\ncontains zero"))+
              
              coord_flip() +
              theme_bw() + 
              theme(axis.text.y = element_text( size=7), 
                    axis.text.x=element_text(size=7),
                    axis.title = element_text(size=7), 
                    legend.text = element_text(size=7)) +
              theme(legend.position=c(0.85, 0.75), 
                    legend.title=element_blank()) +
              xlab(NULL) +
              ylab("Effect on probability\nof receiving seeding")+
              guides(linetype=FALSE, color=FALSE, shape=FALSE) 

postplot2

legend <- get_legend(postplot)


write.csv(posterior, "sourcefigs/fig2b.csv")


##c) Marginal effects plots for covariates from propensity/treatment probability model #################

## visualize differences in distributions of covariates:
dens_road <- ggplot(sagesites, aes(x=distroad, 
                                   color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) +
  scale_fill_manual(values=c("#CC3333", "#197bdd")) +
  xlab("Distance from road") +mytheme

dens_fire <- ggplot(sagesites, aes(x=Shape_Area, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) +
  scale_fill_manual(values=c("#CC3333", "#197bdd")) +
  xlab("Fire size (Ha)") + mytheme

dens_ppt <- ggplot(sagesites, aes(x=wintersprppt, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) +
  scale_fill_manual(values=c("#CC3333", "#197bdd")) + xlim(0,450)+
  xlab("Feb-Apr mean total precip (mm, 1985-2014)") + mytheme 

dens_temp <- ggplot(sagesites, aes(x=febaprtmean.tmean, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) +
  scale_fill_manual(values=c("#CC3333", "#197bdd")) +
  xlab("Feb-Apr mean daily temperature (degrees C, 1985-2014)") + mytheme + 
  theme(legend.position=c(0.25, 0.75))


dens_ARTRpre <- ggplot(sagesites, aes(x=covTSF1pf, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd"),
                     "") +
  scale_fill_manual(values=c("#CC3333", "#197bdd"),
                    "") +
  xlab("Estimated pre-fire sagebrush % cover") + mytheme + theme(legend.position=c(0.75, 0.75),
                                                                 legend.text=element_text(size=7))

dens_ARTRpost <- ggplot(sagesites, aes(x=covTSF0, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) + 
  scale_fill_manual(values=c("#CC3333", "#197bdd")) + 
  xlab("Post-fire surviving Sagebrush % cover") + mytheme


dens_clay <- ggplot(sagesites, aes(x=clay, color=Treatment_status)) +
  geom_density(aes(fill=Treatment_status), alpha=0.4) + 
  scale_color_manual(values=c("#CC3333", "#197bdd")) + 
  scale_fill_manual(values=c("#CC3333", "#197bdd")) + 
  xlab("Soil percent clay (%)") + mytheme



## Marginal effect of precipitation:

n <-nrow(sagesites)

ppt.gradient <- seq(from=min(sagesites$wintersprppt, na.rm=T), 
                    to=max(sagesites$wintersprppt),
                    length.out=100)

preds <- sagesites %>%
  data_grid(wintersprppt = ppt.gradient,
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = 
              mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions = "Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore)

avgprediction <-  preds %>% group_by(wintersprppt) %>%
  median_qi(.value)
avgprediction[which.max(avgprediction$.value),]

write.csv(preds, "sourcefigs/fig3a.csv")

cf.ppt <-   ggplot(preds, aes(x = wintersprppt, y = .value)) +
  stat_lineribbon(
    .width = c(.5, 0.95), alpha = 0.35, 
    fill="#197bdd", color="#0033CC", size=1) +
  ylab("Probability of\nreceiving seeding")+ 
  xlab("Nov-Apr total precipitation (mm)") +
  ylim(0, 1) +
  mytheme + theme(axis.text.x=element_text(size=8))

cf.ppt 


## Marginal effect of temperature
temp.gradient <- seq(from=min(sagesites$febaprtmean.tmean, 
                              na.rm=T), 
                     to= max(sagesites$febaprtmean.tmean, 
                             na.rm=T),
                     length.out=100)

preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = temp.gradient,
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            heatload = mean(heatload, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions ="Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

avgprediction <-  preds %>% group_by(febaprtmean.tmean) %>%
  median_qi(.value)
avgprediction[which.max(avgprediction$.value),]

write.csv(preds, "sourcefigs/fig3b.csv")

cf.temp <-   ggplot(preds, aes(x = febaprtmean.tmean, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), alpha = 0.35, fill="#197bdd", 
                  color="#0033CC", size=1) + 
  ylab("Probability of\nreceiving seeding")+ 
  xlab(expression("Feb-Apr mean temperature " ( degree*C), " ")) +
  ylim(0, 1) +
  mytheme+ theme(axis.text.x=element_text(size=8))

cf.temp



## Marginal effect of prefire cover
cov.gradient <- seq(from=min(1, na.rm=T), 
                    to= max(sagesites$covTSF1pf, na.rm=T),
                    length.out=100)


preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = cov.gradient,
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            heatload = mean(heatload, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions ="Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

avgprediction <-  preds %>% group_by(covTSF1pf) %>%
  median_qi(.value)
avgprediction[which.max(avgprediction$.value),]

write.csv(preds, "sourcefigs/fig3c.csv")

cf.pfcover <-   ggplot(preds, aes(x = covTSF1pf, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), alpha = 0.35, 
                  fill="#197bdd", color="#0033CC", size=1) + 
  ylab("Probability of\nreceiving seeding")+ 
  xlab("Pre-fire sagebrush cover (%)") +
  ylim(0, 1) +
  mytheme

cf.pfcover



## Marginal effect of postfire sagebrush cover
cov.gradient2 <- seq(from=min(sagesites$covTSF0, na.rm=T), 
                     to= max(sagesites$covTSF0, na.rm=T),
                     length.out=100)


preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = cov.gradient2,
            febaprtmean.tmean = mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            heatload = mean(heatload, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            Level3Ecoregions = "Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

write.csv(preds, "sourcefigs/fig3d.csv")

cf.postcover <-   ggplot(preds, aes(x = covTSF0, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), 
                  alpha = 0.35, fill="#197bdd", 
                  color="#0033CC", size=1) + 
  ylab("Probability of\nreceiving seeding")+ 
  xlab("Surviving, unburned sagebrush cover (%)") +
  ylim(0, 1) + 
  mytheme + theme(axis.title.x=element_text(size=7))

cf.postcover


avgprediction <-  preds %>% group_by(covTSF0) %>%
  median_qi(.value)
avgprediction[which.max(avgprediction$.value),]



## Marginal effect of distance from the road
road.gradient <- seq(from=min(sagesites$distroad, na.rm=T), 
                     to= max(sagesites$distroad, na.rm=T),
                     length.out=100)


preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            distroad = road.gradient,
            heatload = mean(heatload, na.rm=T),
            Level3Ecoregions="Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

write.csv(preds, "sourcefigs/fig3e.csv")

cf.road <-   ggplot(preds, aes(x = distroad*28.29, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), alpha = 0.35, 
                  fill="#197bdd", color="#0033CC", size=1) + 
  ylab("Probability of\nreceiving seeding")+ 
  xlab("Distance from a major road (km)") +
  ylim(0, 1) +
  mytheme

cf.road


### Marginal effect of fire size
fire.gradient <- seq(from=min(sagesites$Shape_Area, na.rm=T), 
                     to= max(sagesites$Shape_Area, na.rm=T),
                     length.out=100)


preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = fire.gradient,
            clay = mean(clay, na.rm=T),
            sand = mean(sand, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            heatload = mean(heatload, na.rm=T),
            Level3Ecoregions = "Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

write.csv(preds, "sourcefigs/fig3f.csv")

cf.fire <-   ggplot(preds, aes(x = Shape_Area*0.00012, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), 
                  alpha = 0.35, 
                  fill="#197bdd", 
                  color="#0033CC", size=1) + 
  scale_x_continuous(breaks = seq(0, 100000, by = 50000)) +
  ylim(0, 1) +
  ylab("Probability of\nreceiving seeding")+ xlab("Fire size (hectares)") +
  mytheme

cf.fire



### Marginal effect of percent clay
clay.gradient <- seq(from=min(sagesites$clay, na.rm=T), 
                     to= max(sagesites$clay, na.rm=T),
                     length.out=100)


preds <- sagesites %>%
  data_grid(wintersprppt = mean(wintersprppt, na.rm=T),
            covTSF1pf = mean(covTSF1pf, na.rm=T),
            covTSF0 = mean(covTSF0, na.rm=T),
            febaprtmean.tmean = mean(febaprtmean.tmean, na.rm=T),
            elevation = mean(elevation, na.rm=T),
            Shape_Area = mean(Shape_Area, na.rm=T),
            clay = clay.gradient,
            sand = mean(sand, na.rm=T),
            distroad = mean(distroad, na.rm=T),
            heatload = mean(heatload, na.rm=T),
            Level3Ecoregions = "Snake River\nPlain") %>%
  add_fitted_draws(probtreatbefore, draws=500)

write.csv(preds, "sourcefigs/fig3-clay.csv")

cf.clay <-   ggplot(preds, aes(x = clay, y = .value)) +
  stat_lineribbon(.width = c(.5, 0.95), alpha = 0.35, 
                  fill="#197bdd", color="#0033CC", size=1) + 
  ylab("Probability of\nreceiving seeding")+ xlab("Soil percent clay (%)") +
  ylim(0, 1)+
  mytheme

cf.clay

## ADD density plots to marginal effects plots ###

theme_bare <- theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  #axis.ticks.length = unit(0, "lines"), # Error 
  axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
  panel.background = element_rect(fill = "white"), 
  panel.border = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.margin = unit(c(0,0,0,0), "lines"), 
  plot.background = element_rect(fill = "white"),
  plot.margin = unit(c(0,0,0,0), "lines")
)

legend <- get_legend(dens_ARTRpre)

ppt2 <- cf.ppt %>%
  insert_xaxis_grob(dens_ppt + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw()                  

temp2 <- cf.temp %>%
  insert_xaxis_grob(dens_temp + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 

pfcover2 <- cf.pfcover %>%
  insert_xaxis_grob(dens_ARTRpre + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 

postcover2 <- cf.postcover  %>%
  insert_xaxis_grob(dens_ARTRpost + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 

roads2 <- cf.road %>%
  insert_xaxis_grob(dens_road + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 

fires2 <- cf.fire  %>%
  insert_xaxis_grob(dens_fire + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 

clay2 <- cf.clay  %>%
  insert_xaxis_grob(dens_clay + theme_bare, 
                    grid::unit(0.75, "in"), position = "top") %>%
  ggdraw() 


## Grid plots for manuscript:
cfs.propensity <- plot_grid(cf.pfcover, cf.postcover, 
                            cf.ppt, cf.temp,  
                            cf.clay,  cf.road, 
                            nrow=3, ncol=2, 
                            labels=c("a", "b", "c", "d", "e", "f"),
                            scale = 0.90)

cfs.propensity


parplots <- plot_grid(postplot, postplot2, 
                      nrow=1, ncol=2,  rel_widths = c(1, 0.7))
parplots                      

matchingplot <- plot_grid(prop.grid, parplots, 
                          ncol=1, nrow=2, rel_heights=c(0.5, 1),
                          labels=c("a", "b"))
matchingplot



cfs.dens <- (pfcover2 + labs(tag="a") | postcover2 + labs(tag="b") + inset_element(legend, 0.6, 0.8, 0.9, 0.9)) / 
  (ppt2 + labs(tag="c") | temp2 + labs(tag="d")) / 
  (roads2 + labs(tag="e") | clay2 + labs(tag="f"))


### 5) Save plots ######

pdf('figures/effect of matching_final.pdf',
     width=7.08661, height=7
)
matchingplot
dev.off()


pdf('figures/effect of matching_other.pdf',
     width=7, height=7
)
matchingplot2
dev.off()


pdf('figures/sources of bias_final.pdf',
     width=6.5, height=8
)
cfs.propensity
dev.off()


pdf('figures/sources of bias_density_final.pdf',
     width=5, height=7
)
cfs.dens
dev.off()

rm(list= ls()[!(ls() %in% c('proptreatafter','probtreatbefore', 
                            'summ.stats', 'psmmodel',
                            'matches', 'sagesites', 'pixels'))])


### 6) Assess pre-treatment trends in treated v. untreated pixels ####
timeseries <- select(sagesites, covTSF10pf:covTSF1pf, treatment_status) %>% na.omit()
colnames(timeseries) <- c(-10:-1, "treatment")
timeseries <- gather(timeseries, time, sagecov, "-10":"-1") 


ggplot(timeseries, aes(x=as.numeric(time), y=sagecov, color=as.factor(treatment))) + geom_smooth(method="loess") + 
  scale_color_manual(values=c("#CC3333", "#197bdd"), 
                     labels=c("Untreated Post-fire", "Treated Post-fire")) +
  scale_fill_manual(values=c("#CC3333", "#197bdd"),
                    labels=c("Untreated", "Treated")) +
  xlab("Years before fire")+ ylab("Sagebrush cover") +
  theme_bw() + theme(legend.title=element_blank()) 


