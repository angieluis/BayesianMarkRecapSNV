#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
####### Data format: 
### Capture Histories: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)

###### Covariate dataframes called temporal.covariates & individual.covariates
# temporal.covariates are weekly so that we can simulate z on 
# a weekly time scale. Prim is on the primary occasions
# Use the same (monthly) temporal covariate data for the whole month
# individual covariates include web, sex

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
load("ZunipemaCH.RData")
#load("Zunitemporalcovariates.RData")
#load("Zuniindividualcovariates.RData")

CH.secondary <- Zuni.pema.Ch.secondary
CH.primary <- primary.ch.fun(CH.secondary)
time.int <- Zuni.primary.time.int.weeks 
dates <- Zuni.dates # dates trapped
individual.covariates <- Zuni.individual.covariates


### We don't have temporal covariate data for web 3 so git rid of those data for now.
del <- which(individual.covariates$web==3)
CH.secondary <- lapply(CH.secondary,function(x){x[-del,]})
CH.primary <- CH.primary[-del,]
individual.covariates <- individual.covariates[-del,]
# renumber ID column
individual.covariates$ID <- 1:length(individual.covariates$ID)

temporal.covariates <- weekly.temporaldata.fun(
  dates=dates, # dates trapped
  data= sw.temp.data, # all the temporal covariate data in long format
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(ndvi=0,ndvi=1,tmax=3,tmin=5),# this means I want ndvi for each grid at no lag and t-1 lag, as well as tmax at t-3, and tmin at t-5
  individual.covariates=individual.covariates)
# output is a list with the first element called weekly.longdata to be used for long capture histories below
# other elements of the list (if individual.covariates are included) are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_1 is a matrix of NDVI lag 1 with dimension [individual, month], so can be used for phi~ndvi_1[i,m]


### See model code in "robust_CJS_weekly_phi_NDVI_p_dot_c_dot.bug"


################### Do the data manipulation in R

date()
obs.dat.full <- weekly.longdataCH.fun(CH.secondary, temporal.covariates$weekly.longdata, p_or_c=TRUE)
date() # 22 min

first.caught <- apply(CH.primary,1,function(x){min(which(x>0))})
individual.covariates$f.week <- temporal.covariates$weekly.longdata$week[match(first.caught,temporal.covariates$weekly.longdata$Prim)]



# this creates a weekly capture history to pass into the
# initial values and known state functions
weeklyCH <- weekly.primaryCH.fun(CH.primary,temporal.covariates$weekly.longdata)

##### Bundle data
bugs.data <- list(
  y = obs.dat.full$State,
  prim = obs.dat.full$Prim,
  sec = obs.dat.full$Sec,
  id = obs.dat.full$ID,
  f = individual.covariates$f.week, 
  p.or.c = obs.dat.full$p.or.c,
  web = individual.covariates$web,
  sex = individual.covariates$sex,
  nind = dplyr::n_distinct(obs.dat.full$ID), 
  n.weeks = max(obs.dat.full$week), 
  #time.int = time.int,
  n.obs = nrow(obs.dat.full),
  weeklyCH = weeklyCH,
  z = known.state.cjs(weeklyCH), 
  cjs.init.z=cjs.init.z, 
  CH.primary=CH.primary,
  covariate.month = temporal.covariates$weekly.longdata$month,
  week = obs.dat.full$week,
  ndvi_0 = temporal.covariates$ndvi_0,
  ndvi_1 = temporal.covariates$ndvi_1,
  tmax_3 = temporal.covariates$tmax_3,
  tmin_5 = temporal.covariates$tmin_5
) 

#initial values
inits=function(){list(z=cjs.init.z(weeklyCH,f), alpha.month=runif(12,0,1), mean.p=runif(1,0,1), mean.c=runif(1,0,1), alpha.0=runif(1,0,1), alpha.ndvi_0=runif(1,0,1), alpha.ndvi_1=runif(1,0,1),alpha.tmax_3=runif(1,0,1), alpha.tmin_5=runif(1,0,1) )} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.month","alpha.ndvi_0", "alpha.ndvi_1","alpha.tmax_3","alpha.tmin_5")


date()
Z12.weekly.rcjs.phi.ndvi_1.p.c.constant=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_weekly_phi_ndvi_0_ndvi_1_tmax_3_tmin_5_p_dot_c_dot.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() # 23 hours

save.image("Z12weeklymodels.RData")

Z12.weekly.rcjs.phi.ndvi_1.p.c.constant
hist.plot.fun(Z12.weekly.rcjs.phi.ndvi_1.p.c.constant)
chain.plot.fun(Z12.weekly.rcjs.phi.ndvi_1.p.c.constant)


###############################################################################
### Model Comparisons
## need to use new WAIC instead of DIC
###############################################################################
mod.names <- objects()[grep("Z2.weekly.rcjs",objects())]
mod.list <- list()
mod.table <- data.frame(model=mod.names,npar=rep(NA,length(mod.names)),
                        DIC=rep(NA,length(mod.names)))
for(i in 1:length(mod.names)){
  mod.list[[i]]<- get(mod.names[i]) 
  mod.table$DIC[i] <- mod.list[[i]]$BUGSoutput$DIC
  mod.table$npar[i] <- dim(mod.list[[i]]$BUGSoutput$summary)[1]-1
}
mod.table <- mod.table[order(mod.table$DIC),]
mod.table$delta.DIC <- mod.table$DIC-min(mod.table$DIC)
