#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
####### Data format: 
### Capture Histories: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)


#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
load("Zuni12pemaCH.RData")
load("UNMdata.RData")
sw.temp.data <- read.csv("sw_merge.csv")

CH.secondary <- Zuni12.pema.Ch.secondary
CH.primary <- primary.ch.fun(CH.secondary)
time.int <- Zuni12.primary.time.int.weeks 


#rename secondary occasions from 1 to number of days #prob don't need to do this for Emily's capture history code
for(i in 1:length(CH.secondary)){
  colnames(CH.secondary[[i]]) <- 1:dim(CH.secondary[[i]])[2]
}


# number of secondary occasions
n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))




covariate.data <- monthly.covariate.function(
  capture.data=Zuni12.pema.data, #just pema so no overlap w other sites/species
  CH.secondary=CH.secondary, # list
  tags=rownames(CH.secondary[[1]]),
  sessions=Zuni.sessions, #sessions trapped even if no pema caught
  temporal.data= sw.temp.data, # all the temporal covariate data in long format
  longdata=TRUE,
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(ndvi=0,ndvi=1,tmax=3,tmin=5)# this means I want ndvi for each grid at no lag and t-1 lag, as well as tmax at t-3, and tmin at t-5
)
# output is a list with the first element individual.covariates data frame
# second element called temporal.covariates to be used for long capture histories below
# other elements of the list are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_1 is a matrix of NDVI lag 1 with dimension [individual, month] (where month is longmonth not just months trapped), so can be used for phi~ndvi_1[i,m]. 





## array of dim [indiv,total months, max days] #  those not trapped=NA 
# matches up to p[i,m,d]
p.or.c <- p.or.c.array.fun(CH.secondary, covariate.data$temporal.covariates,n.sec.occ,list=TRUE)



### See model code in "robust_CJS_monthlylist_maxcov.bug"

################### Do all the data manipulation in R to create long data to remove loops in the observation code

obs.dat <- monthly.longdataCH.fun(CH.secondary, covariate.data$temporal.covariates, covariate.data$individual.covariates)
  

monthlyCH <- monthly.primaryCH.fun(CH.primary,covariate.data$temporal.covariates)


##### Bundle data
bugs.data <- list(
  nind = dplyr::n_distinct(obs.dat$ID), 
  n.obs = length(obs.dat$State),
  y = obs.dat$State, #observed yes/no
  prim = obs.dat$Prim, #primary occasion of observation 
  longmonth.obs = obs.dat$long.month, # longmonth of observation
  sec = obs.dat$Sec, #secondary occasion of obs
  id = obs.dat$ID, # individual (row of CH)
  f.longmonth = covariate.data$individual.covariates$f.longmonth, # longmonth first caught, for simulating z
  p.or.c = p.or.c, # an array [i,m,d]to fit into p models
  web = covariate.data$individual.covariates$web,
  sex = covariate.data$individual.covariates$sex,#
  n.sec.occ = n.sec.occ, #
  max.secondary.occasions = max(obs.dat$Sec), 
  n.primary.occasions = max(obs.dat$Prim), 
  time.int = time.int,
  monthlyCH = monthlyCH,
  z = known.state.cjs(monthlyCH),
  cjs.init.z=cjs.init.z,
  #CH.primary=CH.primary,
  n.months = dim(covariate.data$temporal.covariates)[1], #140
  covariate.month = covariate.data$temporal.covariates$month, #1:12
  long.month = covariate.data$temporal.covariates$long.month, #1:140
  months.trapped = covariate.data$temporal.covariates$long.month[which(is.finite(covariate.data$temporal.covariates$Prim))], #
  Prim = covariate.data$temporal.covariates$Prim, #1:127
  ndvi_0 = covariate.data$ndvi_0,
  ndvi_1 = covariate.data$ndvi_1,
  tmax_3 = covariate.data$tmax_3,
  tmin_5 = covariate.data$tmin_5
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
