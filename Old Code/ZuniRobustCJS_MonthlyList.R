#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
####### Data format: 
### Capture Histories: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)

###### Covariate dataframe called temporal.covariates
# temporal.covariates$long.month is every month of the start of the dataset 
  # to the end starting at 1 
# temporal.covariates$covariate.prim is each primary session. If there are
  # missing months in the dataset then there will be more than one long.month
  # to a covariate.prim (goes from 1 to the number of sessions sampled)
# temporal.covariates$month refers to month of the year (seasonality)

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
load("ZunipemaCH.RData")
load("UNMdata.RData")
sw.temp.data <- read.csv("sw_merge.csv")
#load("Zunitemporalcovariates.RData")
#load("Zuniindividualcovariates.RData")

CH.secondary <- Zuni.pema.Ch.secondary
CH.primary <- primary.ch.fun(CH.secondary)
time.int <- Zuni.primary.time.int.weeks 
temporal.covariates <- Zuni.temporal.covariates
individual.covariates <- Zuni.individual.covariates


### We don't have temporal covariate data for Zuni web 3 so git rid of those data for now.
del <- which(individual.covariates$web==3)
CH.secondary <- lapply(CH.secondary,function(x){x[-del,]})
CH.primary <- CH.primary[-del,]
individual.covariates <- individual.covariates[-del,]
# renumber ID column
individual.covariates$ID <- 1:length(individual.covariates$ID)


# use sex as numeric not factor
individual.covariates$sex <- as.numeric(as.character(individual.covariates$sex))
# some individuals have NA for sex. I will split the difference for estimates by setting NAs = 0.5
individual.covariates$sex[which(is.na(individual.covariates$sex))] <- 0.5


#rename secondary occasions from 1 to number of days #prob don't need to do this for Emily's capture history code
for(i in 1:length(CH.secondary)){
  colnames(CH.secondary[[i]]) <- 1:dim(CH.secondary[[i]])[2]
}
# number of secondary occasions
n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))


temporal.covariates <- monthly.temporaldata.fun(
  capture.data=UNMdata,
  temporal.data= sw.temp.data, # all the temporal covariate data in long format
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(ndvi=0,ndvi=1,tmax=3,tmin=5),# this means I want ndvi for each grid at no lag and t-1 lag, as well as tmax at t-3, and tmin at t-5
  individual.covariates=individual.covariates)
# output is a list with the first element called monthly.longdata to be used for long capture histories below
# other elements of the list (if individual.covariates are included) are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_1 is a matrix of NDVI lag 1 with dimension [individual, month], so can be used for phi~ndvi_1[i,m]

## array of dim [indiv,total months, max days] #  those not trapped=NA 
# matches up to p[i,m,d]
p.or.c <- p.or.c.array.fun(CH.secondary, temporal.covariates,n.sec.occ,list=TRUE)



### See model code in "robust_CJS_monthlylist_maxcov.bug"

################### Do all the data manipulation in R to create long data to remove loops in the observation code

obs.dat <- monthly.longdataCH.fun(CH.secondary, temporal.covariates$monthly.longdata, individual.covariates)
  

monthlyCH <- monthly.primaryCH.fun(CH.primary,temporal.covariates$monthly.longdata)


##### Bundle data
bugs.data <- list(
  y = obs.dat$State,
  prim = obs.dat$Prim,
  sec = obs.dat$Sec,
  id = obs.dat$ID,
  f = first_obs$f.longmonth, 
  p.or.c = p.or.c, # an array to fit into p models
  web = individual.covariates$web,
  sex = as.numeric(as.character(individual.covariates$sex)),
  nind = dplyr::n_distinct(obs.dat$ID), 
  n.sec.occ = n.sec.occ,
  max.secondary.occasions = max(obs.dat$Sec), 
  n.primary.occasions = max(obs.dat$Prim), 
  time.int = time.int,
  n.obs = nrow(obs.dat),
  monthlyCH = monthlyCH,
  z = known.state.cjs(monthlyCH),
  cjs.init.z=cjs.init.z,
  #CH.primary=CH.primary,
  n.months = dim(temporal.covariates$monthly.longdata)[1],
  covariate.month = temporal.covariates$monthly.longdata$month,
  months.trapped = temporal.covariates$monthly.longdata$long.month[which(is.finite(temporal.covariates$monthly.longdata$Prim))],
  #Prim = temporal.covariates$monthly.longdata$Prim,
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
Z12.monthlylist.rcjs.maxcov=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_monthlylist_maxcov.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() 
eptm <- proc.time()
eptm-ptm
# completed, took several days, but not sure exactly how long

save.image("Z12monthlylistmodels.RData")

Z12.monthlylist.rcjs.maxcov
hist.plot.fun(Z12.monthlylist.rcjs.maxcov)
chain.plot.fun(Z12.monthlylist.rcjs.maxcov)
