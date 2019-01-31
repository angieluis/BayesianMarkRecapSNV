#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
####### Data format: 
### Capture Histories: array of capture histories 
# with dimensions [individual, primary, secondary] 
# the third dimension (secondary occasions) is the longest of any of the secondary
# occasions. Then also have a vector which gives the number of secondary occasions
# per primary occasion)


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
#sw.temp.data <- read.csv("sw_merge.csv)

CH.secondary <- Zuni.pema.Ch.secondary.array
n.sec.occ <- Zuni.secondary.occasions 
CH.primary <- apply(CH.secondary,c(1,2),sum)
CH.primary <- replace(CH.primary,CH.primary>1,1)

time.int <- Zuni.primary.time.int.weeks 
dates <- Zuni.dates # dates trapped
individual.covariates <- Zuni.individual.covariates
# use sec as numeric not factor
individual.covariates$sex <- as.numeric(as.character(individual.covariates$sex))


### We don't have temporal covariate data for web 3 so git rid of those data for now.
del <- which(individual.covariates$web==3)
CH.secondary <- CH.secondary[-del,,]
CH.primary <- CH.primary[-del,]
individual.covariates <- individual.covariates[-del,]
# renumber ID column
individual.covariates$ID <- 1:length(individual.covariates$ID)

# number of secondary occasions
#n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))

temporal.covariates <- monthly.temporaldata.fun(
  capture.data=UNMdata,
  temporal.data= sw.temp.data, # all the temporal covariate data in long format
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(ndvi=0,ndvi=1,tmax=3,tmin=5),# this means I want ndvi for each grid at no lag and t-1 lag, as well as tmax at t-3, and tmin at t-5
  individual.covariates=individual.covariates)
# output is a list with the first element called monthly.longdata to be used for long capture histories below
# other elements of the list (if individual.covariates are included) are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_1 is a matrix of NDVI lag 1 with dimension [individual, month], so can be used for phi~ndvi_1[i,m]

## 
p.or.c <- p.or.c.array.fun(CH.secondary, temporal.covariates$monthly.longdata,n.sec.occ,list=FALSE)

# create a vector of first marking
f <- apply(CH.primary, 1, function(x) min(which(x!=0)))

### See model code in "robust_CJS_monthly_maxcov.R"

# some individuals have NA for sex. I will split the difference for estimates but setting NAs = 0.5
individual.covariates$sex[which(is.na(individual.covariates$sex))] <- 0.5

# this creates a monthly capture history to pass into the
# initial values and known state functions
monthlyCH <- monthly.primaryCH.fun(CH.primary,temporal.covariates$monthly.longdata)


##### Bundle data
bugs.data <- list(
  y = CH.secondary,
  f = f, # first trap occasion
  f.longmonth = temporal.covariates$monthly.longdata$long.month[match(f,temporal.covariates$monthly.longdata$Prim)], #longmonth first trapped 
  p.or.c = p.or.c, # now an array to into p models
  cjs.init.z = cjs.init.z, 
  CH.primary = CH.primary,
  monthlyCH = monthlyCH,
  nind = dim(CH.primary)[1],
  n.months = dim(temporal.covariates$monthly.longdata)[1],
  covariate.month = temporal.covariates$monthly.longdata$month,
  months.trapped = temporal.covariates$monthly.longdata$long.month[which(is.finite(temporal.covariates$monthly.longdata$Prim))],
  Prim = temporal.covariates$monthly.longdata$Prim,
  n.sec.occ = n.sec.occ,
  web = individual.covariates$web,
  sex = individual.covariates$sex,
  z = known.state.cjs(monthlyCH), 
  ndvi_0 = temporal.covariates$ndvi_0,
  ndvi_1 = temporal.covariates$ndvi_1,
  tmax_3 = temporal.covariates$tmax_3,
  tmin_5 = temporal.covariates$tmin_5
) 

#initial values
inits=function(){list(z=cjs.init.z(monthlyCH,f.longmonth), mean.phi=runif(1,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1),alpha.0=runif(1,0,1),alpha.month=runif(11,0,1),  alpha.ndvi_0=runif(1,0,1), alpha.ndvi_1=runif(1,0,1),alpha.tmax_3=runif(1,0,1), alpha.tmin_5=runif(1,0,1), alpha.male=runif(1,0,1), alpha.month.ndvi_1=runif(11,0,1),sigma.0=runif(1,0,1), sigma.recap=runif(1,0,1),sigma.male=runif(1,0,1),sigma.month=runif(11,0,1) )} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.month","alpha.ndvi_0", "alpha.ndvi_1","alpha.tmax_3","alpha.tmin_5","alpha.male","alpha.month.ndvi_1","sigma.0","sigma.recap","sigma.male","sigma.month")


sptm <- proc.time()
date()
Z12.monthly.rcjs.maxcov=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_monthly_maxcov.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() 
eptm <- proc.time()
eptm-sptm
# completed, took several days, but not sure exactly how long

save.image("Z12monthlyarraymodels.RData")

Z12.monthly.rcjs.maxcov
hist.plot.fun(Z12.monthly.rcjs.maxcov)
chain.plot.fun(Z12.monthly.rcjs.maxcov)


###############################################################################
### Model Comparisons
## need to use new WAIC instead of DIC
###############################################################################
mod.names <- objects()[grep("Z2.monthly.rcjs",objects())]
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
