#### look at qPCR stan model estimates

rm(list=ls())
library(ggplot2)
library(tidyverse)

load("pcod_qPCR_fitted_20240222.RData")

###### Part 1 - teasing apart variability 

### let me start by plotting the sd of the qPCR data (similar figure 4)  

## qPCR standards
sd.among.pcr.stand <- Output.qpcr$sd.among.pcr.stand 

sd.among.pcr.stand %>%
  ggplot() +
  geom_errorbar(aes(x=pcr.rep,ymin=q.05,ymax=q.95),width=0) +
  geom_point(aes(x=pcr.rep, y=mean.SD),size=3) +
  #ylab("") + 
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")

## PCR samples 
sd.among.pcr.samp <- Output.qpcr$sd.among.pcr.samp 

sd.among.pcr.samp %>%
  ggplot() +
  geom_errorbar(aes(x=pcr.rep,ymin=q.05,ymax=q.95),width=0) +
  geom_point(aes(x=pcr.rep, y=mean.SD),size=3) +
  #ylab("") + 
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")

## Bottles
sd.among.bottles <- Output.qpcr$sd.among.bottles 

sd.among.bottles %>%
  ggplot() +
  geom_errorbar(aes(x=site_name,ymin=q.05,ymax=q.95),width=0) +
  geom_point(aes(x=site_name, y=mean.SD),size=3) +
  #ylab("") + 
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")

## PCR + Bottles 
sd.pcr.plus.bottle <- Output.qpcr$sd.pcr.plus.bottle 

sd.pcr.plus.bottle %>%
  ggplot() +
  geom_errorbar(aes(x=site_name,ymin=q.05.TOT.SD,ymax=q.95.TOT.SD),width=0) +
  geom_point(aes(x=site_name, y=TOT.SD),size=3) +
  #ylab("") + 
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")

##  Among sites 
# sd.among.time.given.site <- Output.qpcr$sd.among.time.given.site
# 
# sd.among.time.given.site %>%
#   ggplot() +
#   geom_errorbar(aes(x=month,ymin=q.05,ymax=q.95),width=0) +
#   geom_point(aes(x=month,y=mean.SD),size=3) +
#   #ylab("") + 
#   theme(axis.text.x = element_text(angle = 45),
#         legend.position = "none")

## okay, now create a single data frame to plot all of these 
sd.qPCR <- sd.among.pcr.stand %>%
  mutate(group = "PCR standards") %>%
  mutate(rep = ifelse(pcr.rep == "among pcr", "mean", "individual")) %>%
  rbind(sd.among.pcr.samp %>%
          mutate(group = "PCR samples") %>% 
          mutate(rep = ifelse(pcr.rep == "among pcr", "mean", "individual"))) %>%
  rbind(sd.among.bottles %>%
          mutate(group = "Bottles") %>% 
          mutate(rep = ifelse(site_name == "Avg among bottles", "mean", "individual")) %>%
          rename(pcr.rep = site_name) %>%
          dplyr::select(!month)) %>%
  rbind(sd.pcr.plus.bottle %>% 
          mutate(group = "PCR + Bottles") %>% 
          mutate(rep = ifelse(site_name == "among", "mean", "individual")) %>%
          rename(pcr.rep = site_name) %>%
          rename(mean.SD = TOT.SD) %>%
          rename(q.05 = q.05.TOT.SD) %>%
          rename(q.95 = q.95.TOT.SD) %>%
          dplyr::select(!month)) #%>%
  # rbind(sd.among.time.given.site %>%
  #         mutate(group = "Sites") %>%
  #         rename(pcr.rep = month) %>%
  #         filter(pcr.rep == "Avg among month")%>%
  #         mutate(rep = ifelse(pcr.rep == "Avg among month", "mean", "individual")))

# designate order for x-axis 
sd.qPCR$group <- as.factor(sd.qPCR$group)
sd.qPCR$rep <- as.factor(sd.qPCR$rep)

sd.qPCR$group <- factor(sd.qPCR$group, levels = c("PCR standards","PCR samples","Bottles","PCR + Bottles")) #,"Sites"))

my_plot <- sd.qPCR %>%
  ggplot() +
  geom_pointrange(aes(x=group, y=mean.SD, ymin=q.05,ymax=q.95, col=rep),size=0.2, position = position_jitter(width = 0.2)) + 
  xlab("") + 
  ylab("standard deviation (log10 DNA copy number)") + 
  theme_bw() #+ 
  #theme(axis.text.x = element_text(angle = 45)) 
my_plot

#ggsave("my_figures/qPCR.sd.plot.png", plot = my_plot, width = 6, height = 4, dpi = 300)
