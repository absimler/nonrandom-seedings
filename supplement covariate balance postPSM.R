### Nonrandom treatment locations and sagebrush seeding success##########
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

########################################################

## Supplementary analysis: Assess covariate balance in matched v. unmatched dataset; Difference in means tests ####



### Test for differences between PSM subsample and full sample:
## Precip variable:
ppt.before <- stan_glm(wintersprppt~treatment, 
                data=sagesites, family=Gamma(link = "log"))
plot(ppt.before)

ppt.after <- stan_glm(wintersprppt~treatment, 
                data=matchdat, family=Gamma(link = "log"))
plot(ppt.after)


## Temp variable:
hist(sagesites.NA$febaprtmean.tmean)
temp.before <- stan_lm(febaprtmean.tmean~treatment, 
            data=sagesites, prior=R2(what="mean", location=0.1))
plot(temp.before)

temp.after <- stan_lm(febaprtmean.tmean~treatment, 
                      data=matchdat, prior=R2(what="mean", location=0.1))
plot(temp.after)


## Road variable:
road.before <- stan_glm((distroad+0.000001)~treatment, 
                        data=sagesites, family=Gamma(link = "log"))
plot(road.before)

road.after <- stan_glm((distroad+0.000001)~treatment, 
                      data=matchdat, family=Gamma(link = "log"))
plot(road.after)

# Prefire sage variable:
pfsage.before <- stan_glm.nb(covTSF1pf~treatment,
                        data=sagesites)
plot(pfsage.before)

pfsage.after <- stan_glm.nb(covTSF1pf~treatment, 
                       data=matchdat)
plot(pfsage.after)

# surviving sage variable:
survsage.before <- stan_glm.nb(covTSF0~treatment, 
                          data=sagesites)
plot(survsage.before)

survsage.after <- stan_glm.nb(covTSF0~treatment, 
                         data=matchdat)
plot(survsage.after)


# surviving sage variable:
fire.before <- stan_glm((Shape_Area)~treatment, 
                            data=sagesites, 
                        family=Gamma(link = "log"), chains=2)
plot(fire.before)

fire.after <- stan_glm((Shape_Area)~treatment, 
                           data=matchdat, 
                          family=Gamma(link = "log"), chains=2)
plot(fire.after)


# clay variable:
sagesites$clayper <- (sagesites$clay/100) + 0.00001
matchdat$clayper <- (matchdat$clay/100) + 0.00001

clay.before <- stan_betareg(clayper~treatment, 
                        data=sagesites, 
                       link = "logit", chains=2)
plot(clay.before)

clay.after <- stan_betareg(clayper~treatment, 
                           data=matchdat, 
                           link = "logit", chains=2)
plot(clay.after)


# sand variable:
sagesites$sandper <- (sagesites$sand/100) + 0.0001
matchdat$sandper <- (matchdat$sand/100) + 0.0001

sand.before <- stan_betareg(sandper~treatment, 
                            data=sagesites, 
                            link = "logit", chains=2)
plot(sand.before)

sand.after <- stan_betareg(sandper~treatment, 
                           data=matchdat, 
                           link = "logit", chains=2)
plot(sand.after)


# elevation variable:
elevation.before <- stan_glm((elevation)~treatment, 
                        data=sagesites, 
                        family=Gamma(link = "log"), chains=2)
plot(elevation.before)

elevation.after <- stan_glm((elevation)~treatment, 
                       data=matchdat, 
                       family=Gamma(link = "log"), chains=2)
plot(elevation.after)


# heatload variable:
heatload.before <- stan_glm((heatload)~treatment, 
                             data=sagesites, 
                             family=Gamma(link = "log"), chains=2)
plot(heatload.before)

heatload.after <- stan_glm((heatload)~treatment, 
                            data=matchdat, 
                            family=Gamma(link = "log"), chains=2)
plot(heatload.after)




### Plots for illustrating covariate balance following PSM: ######

#### a) results of difference in means tests before/after ######

posterior <- bind_rows(
                  mcmc_intervals_data(pfsage.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="pfsage"),
                  mcmc_intervals_data(pfsage.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="pfsage"),
                  
                  mcmc_intervals_data(survsage.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="survsage"),
                  
                  mcmc_intervals_data(survsage.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="survsage"),
                  mcmc_intervals_data(ppt.before, 
                          prob_outer=0.99, prob=0.5) %>%
                  mutate(dataset = "before", variable="win-spr precip"),
                  mcmc_intervals_data(ppt.after, prob_outer=0.99,
                                     prob=0.5) %>%
                    mutate(dataset = "after", variable="win-spr precip"),
                  mcmc_intervals_data(temp.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="spr temp"),
                  mcmc_intervals_data(temp.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="spr temp"),
                  mcmc_intervals_data(fire.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="fire"),
                  mcmc_intervals_data(fire.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="fire"),
                  mcmc_intervals_data(road.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="road"),
                  mcmc_intervals_data(road.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="road"),
                  mcmc_intervals_data(clay.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="clay"),
                  mcmc_intervals_data(clay.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="clay"),
                  mcmc_intervals_data(sand.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="sand"),
                  mcmc_intervals_data(sand.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="sand"),
                  mcmc_intervals_data(elevation.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="elevation"),
                  mcmc_intervals_data(elevation.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="elevation"),
                  mcmc_intervals_data(heatload.before, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "before", variable="heatload"),
                  mcmc_intervals_data(heatload.after, prob_outer=0.99,
                                      prob=0.5) %>%
                    mutate(dataset = "after", variable="heatload"))
posterior <- subset(posterior, parameter=="treatment")
posterior <- posterior %>% mutate(reorder=c(20:1))

posterior$nonzero <- NA
posterior$nonzero[posterior$ll>0 & posterior$hh>0] <- "nonzero"
posterior$nonzero[posterior$ll<0 & posterior$hh<0] <- "nonzero"
posterior$nonzero[is.na(posterior$nonzero)] <- "zero"

posterior$varnames <- c("Pre-fire\nsagebrush cover",
              "Pre-fire\nsagebrush cover",
               "Post-fire\nsagebrush cover", 
              "Post-fire\nsagebrush cover", 
               "Nov-Apr\nprecipitation",
              "Nov-Apr\nprecipitation",
               "Feb-Apr\nmean temperature",
              "Feb-Apr\nmean temperature",
              "Fire size", 
               "Fire size", "Distance from\nmajor road",
              "Distance from\nmajor road",
              "Percent clay", "Percent clay", "Percent sand", 
              "Percent sand", "Elevation", "Elevation", "Heatload", "Heatload")

posterior$dataset[posterior$dataset=="after"] <- "Matched subset"
posterior$dataset[posterior$dataset=="before"] <- "Full (unmatched) dataset"



# posterior parameter estimate plots
postplot.meandiff <-  ggplot(subset(posterior), 
                           aes(x = reorder(varnames, reorder), 
                               color=dataset, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75), size = 3/4) +
  scale_color_manual(name="",
                     values = c("#075A63", "grey50")) +
  scale_shape_manual(values=c(16, 17), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y = element_text( size=14), 
        axis.text.x=element_text(size=14),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.position=c(0.8, 0.8),
        legend.title=element_blank(), 
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 16)) +
  xlab(NULL) +
  ylab("Effect of treatment on covariate value\n(Evidence of difference in means, comparing untreated and treated groups)")+
  guides(linetype=FALSE) +
  ylim(-1.1, 1.1) 

postplot.meandiff


#### b) Plot of distribution of covariates before/after #######
afterdat <- matchdat %>% dplyr::select(covTSF0, covTSF1pf, wintersprppt, 
                              febaprtmean.tmean, distroad, Shape_Area,
                              treatment, clay, sand, elevation, heatload)
beforedat <- sagesites %>% 
              dplyr::select(covTSF0, covTSF1pf, wintersprppt, 
                            febaprtmean.tmean, distroad, Shape_Area,
                            treatment, clay, sand, elevation, heatload)


afterdat$dataset <- "Matched\nsubset"
beforedat$dataset <- "Unmatched\ndataset"

alldat <- rbind(afterdat, beforedat)

alldat$Group[alldat$treatment==1] <- "Treated"
alldat$Group[alldat$treatment==0] <- "Untreated"

mytheme <-  theme_bw() + 
  theme(axis.text.y = element_text(size=13.5), axis.text.x= element_text(size=12),
        axis.title = element_text(size=14), legend.text = element_text(size=16))



median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x, prob=0.25), 
             ymax = quantile(x, prob=0.75))  
}

IQR80 <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x, prob=0.15), # 1st quartile
             ymax = quantile(x, prob=0.85))  # 3rd quartile
}

## TEMPERATURE:
p.temp <-ggplot(alldat,aes(x=as.factor(dataset),y=febaprtmean.tmean,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('February-April Mean Temperature (C)  ') + xlab('') +
  ylim(-7, 8) + mytheme + 
  theme(legend.position="none", legend.box = "horizontal") +
  coord_flip()
p.temp

## PRECIP
p.ppt <-ggplot(alldat,aes(x=as.factor(dataset),y=wintersprppt,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) +  
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Nov-April Total Precipitation (mm)  ') + xlab('') +
  ylim(50, 300) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.ppt

## TEMPERATURE:
p.temp2 <-ggplot(alldat,aes(x=as.factor(dataset),y=novjantmean.tmean,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('November-Jan Mean Temperature (C)  ') + xlab('') +
  ylim(-7, 8) + mytheme + 
  theme(legend.position="none", legend.box = "horizontal") +
  coord_flip()
p.temp2

## PRECIP
p.ppt2 <-ggplot(alldat,aes(x=as.factor(dataset),y=novjanppt.pr,
                          fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) +  
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('November-Jan Total Precipitation (mm)  ') + xlab('') +
  ylim(50, 300) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.ppt2

## Roads
p.road <-ggplot(alldat,aes(x=as.factor(dataset),y=distroad,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Distance from major road (km)') + xlab('') +
  ylim(0, 0.5) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.road


## prefire Sagebrush
p.sagepf <-ggplot(alldat,aes(x=as.factor(dataset),y=covTSF1pf,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Pre-fire sagebrush cover') + xlab('') +
  ylim(0, 30) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.sagepf


## postfire Sagebrush
p.sagepost <-ggplot(alldat,aes(x=as.factor(dataset),y=covTSF0,
                             fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC"), name="") + 
  scale_fill_manual(values=c("#CC3333", "#0066CC"), name="") +
  ylab('Surviving, unburned sagebrush cover') + xlab('') +
  ylim(0, 5) + 
  coord_flip() + mytheme
p.sagepost

legend <- get_legend(p.sagepost)

## fire size
p.fire <-ggplot(alldat,aes(x=as.factor(dataset),y=Shape_Area,
                               fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Fire size (hectares)') + xlab('') +
  ylim(0, 800000000) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.fire

## clay
p.clay <-ggplot(alldat,aes(x=as.factor(dataset),y=clay,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Percent clay (%)') + xlab('') +
  #ylim(0, 800000000) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.clay


## sand
p.sand <-ggplot(alldat,aes(x=as.factor(dataset),y=sand,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Percent sand (%)') + xlab('') +
  #ylim(0, 800000000) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.sand

## elevation
p.elevation <-ggplot(alldat,aes(x=as.factor(dataset),y=elevation,
                           fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Elevation (m)') + xlab('') +
  #ylim(0, 800000000) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.elevation

## heatload
p.heatload <-ggplot(alldat,aes(x=as.factor(dataset),y=heatload,
                                fill=Group, color=Group))+
  geom_flat_violin(aes(fill=Group),position=position_nudge(x=.13,y=0),adjust=1.5,trim=FALSE,alpha=.65,colour=NA)+
  stat_summary(fun.data=median_IQR, position=position_dodge(-0.25), geom="linerange", size=1.5, alpha=0.8) +
  stat_summary(fun.data=IQR80, position=position_dodge(-0.25), geom="linerange", size=1, alpha=0.5) + 
  stat_summary(fun="mean",position=position_dodge(-0.25), geom="point", size=2.5) + 
  scale_color_manual(values=c("#CC3333", "#0066CC")) + 
  scale_fill_manual(values=c("#CC3333", "#0066CC")) +
  ylab('Heatload') + xlab('') +
  #ylim(0, 800000000) + 
  coord_flip() + mytheme + theme(legend.position="none")
p.heatload


beforeafter1 <- plot_grid(p.sagepf, p.sagepost+theme(legend.position = "none"), 
                          p.ppt, p.temp,legend,
                            nrow=3, ncol=2, 
                         labels=c("A", "B", "C", "D"),
                            scale = 0.99)
beforeafter1

beforeafter2 <- plot_grid(p.road, p.fire, p.clay, p.sand, p.elevation,
                          p.heatload, legend,
                          nrow=4, ncol=2, 
                          labels=c("G", "H", "I", "J", "K", "L"),
                          scale = 0.99)
beforeafter2


### SAVE FIGURES ####
jpeg('figures/diffinmeans_beforeafter.jpeg',
     width=9, height=7, units="in", res=600
)
postplot.meandiff
dev.off()

jpeg('figures/matching_beforeafter.jpeg',
     width=9, height=11, units="in", res=600
)
beforeafter1
dev.off()

jpeg('figures/matching_beforeafter2.jpeg',
     width=9, height=11, units="in", res=600
)
beforeafter2
dev.off()



