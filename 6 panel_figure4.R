### Nonrandom treatment locations and sagebrush seeding success #############
##### Simler-Williamson & Germino 2022 #####
#### Contact: allisonsimlerwil@boisestate.edu #####
#### Final submission version, April 2022 ######

### Figure for treatment effect over time, from panel model:
treatmenteff.panel <- readRDS( "modelfits/treatmenteff.panelsitefireid_feb28.rda")

modfit <- posterior_samples(treatmenteff.panel, "^b")

## E{IF treated}
year0.t <- as.data.frame(c(exp(modfit[1] + 
                                 modfit[2] ), "treated", 0))
year1.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[4]), "treated", 1))
year2.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[5]), "treated", 2))
year3.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[6]), "treated", 3))
year4.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[7]), "treated", 4))
year5.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[8]), "treated", 5))
year6.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[9]), "treated", 6))
year7.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[10]), "treated", 7))
year8.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[11]), "treated", 8))
year9.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[12]), "treated", 9))
year10.t <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[3] + modfit[13]), "treated", 10))

treated.preds <- dplyr::bind_rows(year0.t, year1.t, year2.t, year3.t, year4.t, year5.t, year6.t, year7.t, year8.t, year9.t, year10.t) %>%
                  mutate_all( ~replace(., is.na(.), 0)) %>%
                   mutate(year = rowSums(.[3:13])) %>%
                  rename(treatmentdummy = X.treated., sagecov =b_Intercept) %>%
                    dplyr::select(sagecov, treatmentdummy, year)


# E[IF untreated]
year0.ut <- as.data.frame(c(exp(modfit[1] + modfit[2]), "untreated", 0))
year1.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[4]), "untreated", 1))
year2.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[5]), "untreated", 2))
year3.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[6]), "untreated", 3))
year4.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[7]), "untreated", 4))
year5.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[8]), "untreated", 5))
year6.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[9]), "untreated", 6))
year7.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[10]), "untreated", 7))
year8.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[11]), "untreated", 8))
year9.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[12]), "untreated", 9))
year10.ut <- as.data.frame(c(exp(modfit[1] + modfit[2] + modfit[13]), "untreated", 10))

untreated.preds <- dplyr::bind_rows(year0.ut, year1.ut, year2.ut, year3.ut, year4.ut, year5.ut, 
                                  year6.ut, year7.ut, year8.ut, year9.ut, year10.ut) %>%
  mutate_all( ~replace(., is.na(.), 0)) %>%
  mutate(year = rowSums(.[3:13])) %>%
  rename(treatmentdummy = X.untreated., sagecov = b_Intercept) %>%
  dplyr::select(sagecov,treatmentdummy, year)

preds <- rbind(untreated.preds, treated.preds)


## Summarize yhats for each time point

means <- preds %>%
  group_by(treatmentdummy, year) %>%
  mean_qi(.width=c(0.99))


paneltime <-  ggplot2::ggplot(means, ggplot2::aes(x=as.integer(year), y=sagecov,
                                    color=treatmentdummy)) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=.lower, 
                    ymax=.upper), width=.5, size=0.8) +
  
  geom_line(aes(x=as.integer(year), y=sagecov), size=0.8, alpha=0.5, linetype="dashed")+
  scale_color_manual(values=c("#CC3333", "#197bdd"),
                     name="", labels=c("E[If Treated]", "E[If Untreated]")) +
  scale_x_continuous(breaks=seq(0,10,by=1))+
  scale_y_continuous(breaks=seq(0,10,by=1))+
  ylab("Expected sagebrush cover at treated sites")+ 
  xlab("Years since treatment")+ 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position=c(0.15, 0.85),
        axis.text.y = element_text( size=8), 
              axis.text.x= element_text(size=8),
              axis.title = element_text(size=8),
        legend.text = element_text(size=8))

paneltime

ggplot2::ggplot(preds, ggplot2::aes(x=year, y=sagecov,
                                    color=treatmentdummy, fill=treatmentdummy)) +
  stat_lineribbon(.width=c(0.95)) +
  scale_color_manual(values=c("#CC3333", "#197bdd"),
                     name="") +
  scale_fill_manual(values=c("grey60", "lightblue"),
                    name="")+
  #geom_boxplot()+
  ylab("Expected sagebrush cover at treated sites")+ 
  xlab("Years since treatment")+
  theme_bw()


jpeg('figures/treatmenteffect_paneltime.jpeg',
     width=8, height=5, units="in", res=600
)
paneltime
dev.off()

