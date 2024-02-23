### Fitting qPCR statistical model
rm(list=ls())
library(dplyr)
library(mvtnorm)
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
select <- dplyr::select

base.dir <- "~/gadid_aquaria/qPCR_pollock" # INSERT WORKING DIRECTORY HERE.
plot.dir <- "~/gadid_aquaria/qPCR_pollock/my_figures/" # INSERT DIRECTORY TO PLACE PLOT HERE
setwd(base.dir)

QPCR.DAT <- readRDS("./data/qPCR_data_pollock_20240222.RDS")
QPCR.DAT.OLD <- readRDS("./data/qPCR_data_pollock_20240216.RDS")

dat.pcr.control <- QPCR.DAT[[3]]  # Data frame for pcr controls
dat.stand <- QPCR.DAT[[1]]       # Data frame for GACH of known dna concentration (standards)
dat.samp <- QPCR.DAT[[2]]       # qPCR results for field samples.

dat.pcr.control.old <- QPCR.DAT.OLD[[3]]  # Data frame for pcr controls
dat.stand.old <- QPCR.DAT.OLD[[1]]       # Data frame for GACH of known dna concentration (standards)
dat.samp.old <- QPCR.DAT.OLD[[2]]       # qPCR results for field samples.


#### need to do some modifications to dat.samp 
dat.samp <- dat.samp %>%
  mutate(month = 1) %>%         ## add month column to match framework 
  rename(lab_label = lab_lab)
  
dat.stand$density <- as.numeric(dat.stand$density)


#################################################################
#################################################################
# Construct the indices and helper files to make running the statistical model easy.
#################################################################
#################################################################
SITES  <- data.frame(site_name=sort(unique(dat.samp$site_name)),site_idx=1:length(unique(dat.samp$site_name)))
N_site <- max(SITES$site_idx)
MONTH  <- data.frame(month=sort(unique(dat.samp$month)),month_idx=1:length(unique(dat.samp$month)))
N_month <- nrow(MONTH)
BOTTLES <- data.frame(lab_label=sort(unique(dat.samp$lab_label)),bottle_idx=1:length(unique(dat.samp$lab_label)))
N_bottle <- nrow(BOTTLES)

bottles_labels <- dat.samp %>% dplyr::select(month,lab_label,site_name) %>% 
                      group_by(lab_label,month,site_name) %>% summarise(x=length(lab_label)) %>%
                      as.data.frame()

BOTTLES <- left_join(BOTTLES, bottles_labels %>% dplyr::select(month,lab_label,site_name))

PCR     <- data.frame(qpcr_date=sort(unique(dat.samp$qpcr_date)),pcr_idx=1:length(unique(dat.samp$qpcr_date)))
N_pcr   <- nrow(PCR)

SITE_MONTH <- dat.samp %>% group_by(site_name,month) %>% summarise(n_obs = length(site_name)) %>% dplyr::select(-n_obs) %>% as.data.frame()
SITE_MONTH$site_month_idx <- 1:nrow(SITE_MONTH)

SITE_MONTH_BOTTLES <- full_join(BOTTLES,SITE_MONTH)

SITE_MONTH <- full_join(SITE_MONTH,SITES,by="site_name")
SITE_MONTH <- full_join(SITE_MONTH,MONTH,by="month")
N_site_month <- nrow(SITE_MONTH)

### Combine the indices
dat.samp <- left_join(dat.samp,SITES,by="site_name") 
dat.samp <- left_join(dat.samp,MONTH,by="month") 
dat.samp <- left_join(dat.samp,BOTTLES,by="lab_label") 
dat.samp <- left_join(dat.samp,PCR,by="qpcr_date")

dat.samp <- dat.samp %>%
  select(!site_name.y) %>%
  select(!month.y) %>%
  rename(site_name = site_name.x) %>%
  rename(month = month.x)

dat.samp <-  left_join(dat.samp,SITE_MONTH_BOTTLES,by=c("site_name","month","lab_label"))

dat.samp <- dat.samp %>%
  select(!bottle_idx.y) %>%
  rename(bottle_idx = bottle_idx.x)

##### MAKE FILES THAT ARE EASY TO ROLL INTO STAN.
dat.samp.bin <- dat.samp %>% filter(species=="POLLOCK")
dat.samp.count <- dat.samp %>% filter(pres==1, species=="POLLOCK")  ###################################################

dat.stand.bin <- left_join(dat.stand,PCR,by="qpcr_date")
dat.stand.count <- dat.stand.bin %>% filter(pres==1)

SPECIES <- "POLLOCK"
if(SPECIES == "POLLOCK"){
  N_obs_bin   <- nrow(dat.stand.bin)
  N_obs_count <- nrow(dat.stand.count)
}

##### How many bottles and sites-month combinations have >0 counts for the qPCR?
n.zeros.bottle <- dat.samp.bin %>% group_by(site_name,month,lab_label,site_month_idx) %>% 
            summarize(n.tot = length(pres),n.obs = sum(pres)) %>% mutate(frac = n.obs/n.tot) %>% as.data.frame()
n.zeros.site.month <- n.zeros.bottle %>% group_by(site_name,month,site_month_idx) %>% 
            summarize(N.bot = length(lab_label),N.tot = sum(n.tot),N.pos=sum(n.obs)) %>% as.data.frame()
length(which(n.zeros.site.month$N.obs==0))

counter <- SITE_MONTH %>% group_by(month_idx) %>% summarize(counter= length(month)) %>% arrange(month_idx) %>% as.data.frame()
##################################################################
#### MAKE DATA FOR STAN
##################################################################

#OFFSET = -4.5 # Value to improve Fitting in STAN
OFFSET = 0

stan_data = list(

    "bin_stand"   = dat.stand.bin$pres,
    "count_stand" = dat.stand.count$Ct,
    "D_bin_stand" = log10(dat.stand.bin$density),
    "D_count_stand" = log10(dat.stand.count$density),
    
    "bin_samp"    = dat.samp.bin$pres,
    "count_samp"  = as.numeric(dat.samp.count$Ct),

    # Indices and counters
    "N_site"   = N_site,   # Number of Sites
    "N_month"  = N_month,  # Number of months observed
    "N_bottle" = N_bottle, # Number of individual bottles observed.
    "N_pcr"    = N_pcr,    # Number of PCR plates
    "N_site_month" = N_site_month,
    "N_bin_stand"   = length(dat.stand.bin$pres),
    "N_count_stand" = length(dat.stand.count$Ct),
    "N_bin_samp"   = length(dat.samp.bin$pres),
    "N_count_samp" = length(dat.samp.count$Ct),
    
    # Indices for Standards
    "pcr_stand_bin_idx"   = dat.stand.bin$pcr_idx,
    "pcr_stand_count_idx" = dat.stand.count$pcr_idx,
    "pcr_samp_bin_idx"   = dat.samp.bin$pcr_idx,
    "pcr_samp_count_idx" = dat.samp.count$pcr_idx,
    
    # Indices for site-months and bottles
    "site_month_idx" = SITE_MONTH_BOTTLES$site_month_idx,
    "bottle_idx"     = SITE_MONTH_BOTTLES$bottle_idx,
    
    # Index used in calculating Monthly abundance index over space
    "gamma_idx" = SITE_MONTH$month_idx,
    #"counter" =  counter$counter,
    
    # Indices for Samples
    "site_bin_idx"   = dat.samp.bin$site_idx,
    "site_count_idx" = dat.samp.count$site_idx,
    "month_bin_idx"    = dat.samp.bin$month_idx,
    "month_count_idx"  = dat.samp.count$month_idx,
    "bottle_bin_idx"   = dat.samp.bin$bottle_idx,
    "bottle_count_idx" = dat.samp.count$bottle_idx,
    "site_month_bin_idx"   = dat.samp.bin$site_month_idx,
    "site_month_count_idx" = dat.samp.count$site_month_idx,
    
    #Offset of density for improving fitting characteristics
    "OFFSET" = OFFSET
)


 stan_pars = c(
   "beta_0", # intercept for standards
   "beta_1", # slope for standards
   "phi_0",  # logit intercept for standards
   "phi_1",  # logit slope for standard,

   "gamma",  # site-month combinations for log-densities
   "delta",  # random effect for each bottle around the site-month mean.
   
   "D",     # Latent variable for Log-Density in each bottle-site-time
   
   "sigma_stand_int", # variability among standards regression.
   "sigma_pcr",     # variability among samples, given individual bottle, site, and month 
   "tau_bottle"   # variability among bottles, given site, and month
   
)   
    
### INTIAL VALUES

 stan_init_f1 <- function(n.chain,N_bottle,N_pcr,N_site_month){ 
   A <- list()
   for(i in 1:n.chain){
     A[[i]] <- list(
       
       sigma_stand_int = runif(1,0.01,2),
       beta_0 = runif(N_pcr,20,30),
       beta_1 = rnorm(N_pcr,-3,1),
       
       phi_0  = runif(N_pcr,0,20),
       phi_1  = rnorm(N_pcr,5,1),
       D      = rnorm(N_site_month,-4,1),
       gamma  = rnorm(N_site_month,-4,1),
       sigma_pcr = runif(1,0.01,0.4),
       tau_bottle = runif(1,0.01,0.2)
     )
   }
   return(A)
 }
 
 #################################################################### 
 #################################################################### 
 ##### STAN
 #################################################################### 
 #################################################################### 
 N_CHAIN = 5
 Warm = 5000
 Iter = 10000
 Treedepth = 11
 Adapt_delta = 0.90
 

 stanMod = stan(file = 'qPCR_Dryad_Shelton.stan',data = stan_data, 
                verbose = FALSE, chains = N_CHAIN, thin = 5, 
                warmup = Warm, iter = Warm + Iter, 
                control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                pars = stan_pars,
                boost_lib = NULL,
                init = stan_init_f1(n.chain=N_CHAIN,
                                     N_bottle=N_bottle,
                                     N_pcr=N_pcr,
                                     N_site_month=N_site_month
                                    ))
 
 #################################################################### 
 #################################################################### 
 #################################################################### 
    
 pars <- rstan::extract(stanMod, permuted = TRUE)
 # get_adaptation_info(stanMod)
 samp_params <- get_sampler_params(stanMod)
 #samp_params 
 stanMod_summary <- summary(stanMod)$summary
 round(stanMod_summary,2)
 
 base_params <- c(
 "beta_0",
 "beta_1",
 "phi_0",
 "phi_1",
 "sigma_stand_int", # variability among standard regression.
 "sigma_pcr",     # variability among samples, given individual bottle, site, and month 
 "tau_bottle"   # variability among bottles, given site, and month
) 
 
##### MAKE SOME DIAGNOSTIC PLOTS

print(traceplot(stanMod,pars=c("lp__",base_params),inc_warmup=FALSE))
 
#pairs(stanMod, pars = c(base_params), log = FALSE, las = 1)
 
B1 <- apply(pars$beta_0,2,mean)
B2 <- apply(pars$beta_1,2,mean)

P0 <- apply(pars$phi_0,2,mean)
P1 <- apply(pars$phi_1,2,mean)

V0 <-  mean(pars$sigma_stand_int)

# Plot regression against Standard

# X <- seq(-7,0,length.out=1000) - OFFSET
 X <- seq(1, 5, length.out = 1000)
 Y <- t(B1 + B2 %*% t(X ))
 
 STAND.REG <- data.frame(X=X,Y=Y)
 STAND.REG <- melt(STAND.REG,id.vars="X",value.name="Y")
 
 x.lim=c(min(X),max(X))
 y.lim=c(20,40)
 for(i in 1:ncol(Y)){
  plot(Y[,i]~X,xlim=x.lim,ylim=y.lim,type="l",col=2)
  par(new=T)
 }
 plot(dat.stand.count$Ct~log10(dat.stand.count$density),xlim=x.lim,ylim=y.lim)

 
 #BREAKS <- c(-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5)
 BREAKS <- c(1,2,3,4,5)
 LABS   <- BREAKS + OFFSET

lab.temp<- data.frame(variable=c("Y.1","Y.2","Y.3"), 
                      qpcr_date=c("20230213_1","20230213_2","20240213_3"))
STAND.REG <-left_join(STAND.REG,lab.temp)
  

stand.plot <- ggplot(dat.stand) +
   geom_point(aes(x=log(density,10)- OFFSET ,y=Ct,shape=as.factor(qpcr_date)),alpha=0.75) +
   scale_shape_discrete(name="qPCR Date") +
   theme_bw()
stand.plot <- stand.plot +
    geom_line(data=STAND.REG,aes(x=X,y=Y,color=qpcr_date)) +
   scale_color_discrete(name="qPCR Date") +
   ylab("PCR cycle") +
    xlab("DNA copy number") +
    scale_x_continuous(breaks=BREAKS,labels = paste0("1e",LABS),limits = c(0,6))

stand.plot
 
 # Plot occurrence of standard
 Y <- t(P0 + P1 %*% t(X ))
 LOGIT <- data.frame(X=X,Y=plogis(Y))
 LOGIT <- melt(LOGIT,id.vars="X",value.name="Y")
 
 LOGIT <-left_join(LOGIT,lab.temp)
 
stand.plot.pres <- ggplot(dat.stand.bin) +
   geom_jitter(aes(x=log(density,10)-OFFSET,y=pres,shape=as.factor(qpcr_date)),alpha=0.75,width=0,height=0.05) +
   geom_line(data=LOGIT,aes(y=Y,x=X,color=qpcr_date)) +
   theme_bw() +
   scale_shape_discrete(name="qPCR Date") +
   scale_color_discrete(name="qPCR Date") +
   xlab("DNA copy number") +
   ylab("Amplification success") +
   scale_y_continuous(breaks=c(0,1),labels = c("No","Yes"))+
   scale_x_continuous(breaks=BREAKS,labels = paste0("1e",LABS),limits = c(0,6))
 
stand.plot.pres
 
  #setwd(base.dir)
  #setwd(plot.dir)
  #pdf("PCR diagnostics.pdf",width=8,height=7)
  #  print(pcod.stand.plot)  
  #  print(pcod.stand.plot.pres)
  #dev.off()
  
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ##### Extract data of interest, save to file for use elsewhere
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  
  site.month.out <- data.frame(site_month_idx= 1:ncol(pars$gamma), Mean=apply(pars$gamma,2,mean),
                               Sd=apply(pars$gamma,2,sd),
                               data.frame(t(apply(pars$gamma,2,quantile,probs=c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)))))
  bottles.out     <- data.frame(bottle_idx=BOTTLES$bottle_idx, Mean = apply(pars$D,2,mean),
                                Sd=apply(pars$D,2,sd),
                                data.frame(t(apply(pars$D,2,quantile,probs=c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)))))

  SITE_MONTH_summary <- left_join(SITE_MONTH,site.month.out,by="site_month_idx") %>% group_by(site_name,month) %>%
                            summarize(MEAN = mean(Mean),SD = mean(Sd),
                                      q.025= mean(X2.5.),
                                      q.05 = mean(X5.),
                                      q.10 = mean(X10.),
                                      q.25 = mean(X25.),
                                      Median=  mean(X50.),
                                      q.75 = mean(X75.),
                                      q.90 = mean(X90.),
                                      q.95 = mean(X95.),
                                      q.975= mean(X97.5.))
  BOTTLES_summary    <- left_join(BOTTLES,bottles.out, by= "bottle_idx")  %>%
                              rename( q.025= X2.5.,
                                      q.05 = X5.,
                                      q.10 = X10.,
                                      q.25 = X25.,
                                      Median=  X50.,
                                      q.75 = X75.,
                                      q.90 = X90.,
                                      q.95 = X95.,
                                      q.975= X97.5.) 
  
  ####
  sd.of.count <- data.frame(
                        id = c("among.standard","among.pcr"),
                        log.sd.est = c(mean(pars$sigma_stand_int), mean(pars$sigma_pcr)),
                        log.sd.uncert = c(sd(pars$sigma_stand_int), sd(pars$sigma_pcr))
                        )
  
  site_month_samps <- data.frame(site_month_idx=1:N_site_month , t( pars$gamma))
  site_month_samps <- melt(site_month_samps,id.vars="site_month_idx",variable.name="rep") %>% 
    full_join(SITE_MONTH %>% dplyr::select(site_name,month,site_month_idx),.,by="site_month_idx")

### only one time point    
#  sd.among.site.given.time <- rbind(site_month_samps %>% group_by(site_name,rep) %>% summarize(SD=sd(value)) %>% 
#                                      group_by(site_name) %>% summarize(mean.SD=mean(SD),q.05=quantile(SD,probs=0.05),q.95=quantile(SD,probs=0.95)),
#                                  site_month_samps %>% group_by(site_name,rep) %>% summarize(SD=sd(value)) %>% 
#                                      group_by(rep) %>% summarize(val=mean(SD)) %>% mutate(site_name="Avg among site") %>% group_by(site_name) %>%
#                                      summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))
 
  sd.among.time.given.site <- rbind(site_month_samps %>% group_by(month,rep) %>% summarize(SD=sd(value)) %>% 
                                group_by(month) %>% summarize(mean.SD=mean(SD),q.05=quantile(SD,probs=0.05),q.95=quantile(SD,probs=0.95)),
                              site_month_samps %>% group_by(month,rep) %>% summarize(SD=sd(value)) %>% 
                                group_by(rep) %>% summarize(val=mean(SD)) %>% mutate(month="Avg among month") %>% group_by(month) %>%
                                summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))

  bottle_samps <- data.frame(bottle_idx=BOTTLES$bottle_idx, t(pars$D))
  bottle_samps <- melt(bottle_samps,id.vars="bottle_idx",variable.name="rep") %>% 
                    full_join(BOTTLES %>% dplyr::select(site_name,month,bottle_idx),.,by="bottle_idx")
  
  sd.among.bottles <- rbind(bottle_samps %>% group_by(site_name,month,rep) %>% summarize(SD=sd(value)) %>% 
                              group_by(site_name,month) %>% summarize(mean.SD=mean(SD),q.05=quantile(SD,probs=0.05),q.95=quantile(SD,probs=0.95)) %>%
                              mutate(month=as.character(month)),
                            bottle_samps %>% group_by(site_name,month,rep) %>% summarize(SD=sd(value)) %>% 
                              group_by(rep) %>% summarize(val=mean(SD)) %>% mutate(site_name="Avg among bottles",month="all") %>% group_by(site_name,month) %>%
                              summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))

  sd.bottle.model          <- data.frame(mean.SD=mean(pars$tau_bottle),q.05=quantile(pars$tau_bottle,probs=0.05),q.95=quantile(pars$tau_bottle,probs=0.95))
  
  # Using the math of random variables to calculate the expected variability among samples due to sample processing and PCR replicates.
  Y <- seq(30,40,length.out=1000)
  b_0 <- apply(pars$beta_0,2,mean)
  b_1 <- apply(pars$beta_1,2,mean)
  stand.sd <- mean(pars$sigma_stand_int)
  pcr.sd  <- mean(pars$sigma_pcr)
  tot.var <- mean(pars$sigma_stand_int^2 + pars$sigma_pcr^2)
  
  
  ### FIX TO REFLECT SAMPLING UNCERTAINTY.
  X.mean <- (Y - b_0) / b_1
  X.sd.stand <- sqrt(apply(pars$beta_1^(-2),2,"*",pars$sigma_stand_int^2))
  X.sd.pcr   <- sqrt(apply(pars$beta_1^(-2),2,"*",pars$sigma_pcr^2))
  X.sd.tot   <- sqrt(apply(pars$beta_1^(-2),2,"*",pars$sigma_stand_int^2 + pars$sigma_pcr^2))
  
  sd.among.pcr.samp <- rbind(X.sd.pcr %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                                group_by(pcr.rep) %>% summarize(mean.SD = mean(value),q.05=quantile(value,probs=0.05),q.95=quantile(value,probs=0.95)),
                            X.sd.pcr %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                                group_by(id) %>% summarize(val = mean(value)) %>% mutate(pcr.rep="among pcr") %>% group_by(pcr.rep) %>%
                                summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))
  
  sd.among.pcr.stand <- rbind(X.sd.stand %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                            group_by(pcr.rep) %>% summarize(mean.SD = mean(value),q.05=quantile(value,probs=0.05),q.95=quantile(value,probs=0.95)),
                        X.sd.stand %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                            group_by(id) %>% summarize(val = mean(value)) %>% mutate(pcr.rep="among pcr") %>% group_by(pcr.rep) %>%
                            summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))

  sd.among.pcr.tot <- rbind(X.sd.tot %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                              group_by(pcr.rep) %>% summarize(mean.SD = mean(value),q.05=quantile(value,probs=0.05),q.95=quantile(value,probs=0.95)),
                            X.sd.tot %>% as.data.frame() %>% mutate(id=1:nrow(X.sd.pcr)) %>% melt(id.vars="id",variable.name="pcr.rep") %>% 
                              group_by(id) %>% summarize(val = mean(value)) %>% mutate(pcr.rep="among pcr") %>% group_by(pcr.rep) %>%
                              summarize(mean.SD = mean(val),q.05=quantile(val,probs=0.05),q.95=quantile(val,probs=0.95)))

  
  sd.pcr.plus.bottle <- left_join( bottle_samps %>% group_by(site_name,month,rep) %>% summarize(bottle.VAR=var(value)), 
                                X.sd.tot^2 %>% as.data.frame() %>% mutate(rep=paste0("X",1:nrow(X.sd.pcr))) %>% 
                                  melt(id.vars="rep",variable.name="pcr.var") %>%
                                  group_by(rep) %>% summarize(pcr.VAR = mean(value))) %>% 
                                  mutate(tot.var = bottle.VAR + pcr.VAR,tot.sd=sqrt(tot.var)) %>% 
                                  group_by(site_name,month) %>% summarize(TOT.SD = mean(tot.sd), q.05.TOT.SD = quantile(tot.sd,probs=0.05), q.95.TOT.SD = quantile(tot.sd,probs=0.95))
    
  #### add in among site-time average for pcr + bottle error                            
   sd.pcr.plus.bottle <- left_join( bottle_samps %>% group_by(site_name,month,rep) %>% summarize(bottle.VAR=var(value)), 
                  X.sd.tot^2 %>% as.data.frame() %>% mutate(rep=paste0("X",1:nrow(X.sd.pcr))) %>% 
                    melt(id.vars="rep",variable.name="pcr.var") %>%
                    group_by(rep) %>% summarize(pcr.VAR = mean(value))) %>%
                    mutate(tot.var = bottle.VAR + pcr.VAR,tot.sd=sqrt(tot.var)) %>% 
                    group_by(rep) %>% summarize(val = mean(tot.sd)) %>% mutate(site_name="among",month=99) %>%
                    group_by(site_name,month) %>% summarize(TOT.SD = mean(val),  q.05.TOT.SD = quantile(val,probs=0.05), q.95.TOT.SD = quantile(val,probs=0.95)) %>%
                    #dplyr::select(site_name,month,TOT.SD,q.05.TOT.SD,q.95.TOT.SD) %>%
                    rbind(sd.pcr.plus.bottle,.)
  
  #########################################
  # Output from model fits plus basic raw data.  
  Output.qpcr <- list(stanMod = stanMod, stanMod_summary = stanMod_summary,samp = pars, samp_params=samp_params,
                 dat.samp=dat.samp, 
                 dat.samp.bin = dat.samp.bin, 
                 dat.samp.count = dat.samp.count,
                 dat.stand.bin =dat.stand.bin, 
                 dat.stand.count = dat.stand.count,
                 dat.pcr.control.pcod = dat.pcr.control,
                 OFFSET = OFFSET,
                 base_params =base_params,
                 SITE_MONTH = SITE_MONTH, SITE_MONTH_summary = SITE_MONTH_summary,
                 BOTTLES = BOTTLES, BOTTLES_summary = BOTTLES_summary,
                 N_site   = N_site,   # Number of Sites
                 N_month  = N_month,  # Number of months observed
                 N_bottle = N_bottle, # Number of individual bottles observed.
                 N_pcr    = N_pcr,    # Number of PCR plates
                 N_site_month = N_site_month,
                 sd.among.pcr.stand = sd.among.pcr.stand,
                 sd.among.pcr.samp = sd.among.pcr.samp,
                 sd.among.pcr.tot = sd.among.pcr.tot,
                 sd.among.bottles = sd.among.bottles,
                 sd.among.bottles.model = sd.bottle.model,
                 sd.pcr.plus.bottle = sd.pcr.plus.bottle,
                 #sd.among.site.given.time = sd.among.site.given.time,
                 sd.among.time.given.site = sd.among.time.given.site,
                 ### pre-compiled plots
                 stand.plot = stand.plot,
                 stand.plot.pres = stand.plot.pres
                 )
 
  setwd(base.dir)
  save(Output.qpcr,file="pcod_qPCR_fitted_20240222.RData")
