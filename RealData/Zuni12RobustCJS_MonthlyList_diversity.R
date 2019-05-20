#################################################################
# test on smaller dataset

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

load("Zuni12pemaCH.RData")
load("DiversityDataZuni12.RData")
session.list <- Zuni.session.list
CH.secondary <- Zuni12.pema.Ch.secondary
sw.temp.data <- read.csv("sw_merge.csv")
CH.primary <- primary.ch.fun(CH.secondary)


#rename secondary occasions from 1 to number of days #prob don't need to do this for Emily's capture history code
for(i in 1:length(CH.secondary)){
  colnames(CH.secondary[[i]]) <- 1:dim(CH.secondary[[i]])[2]
}


# number of secondary occasions
n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))



# check and make sure using appropriate data
covariate.data <- monthly.covariate.function(
  capture.data=Zuni12.pema.data, #just pema so no overlap w other sites/species
  CH.secondary=CH.secondary, # list
  tags=rownames(CH.secondary[[1]]),
  sessions=session.list$all.sessions, #sessions trapped even if no pema caught
  temporal.data= sw.temp.data, # all the temporal covariate data in long format
  diversity.data= MNAs.diversityZ12,
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(ndvi=2,pt=0, pt=1,ShannonH=0)# this means I want ndvi for each grid at t-2 lag, as well as MNA of P. truei at no lag
)
# output is a list with the first element individual.covariates data frame
# second element called temporal.covariates to be used for long capture histories below
# other elements of the list are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_2 is a matrix of NDVI lag 2 with dimension [individual, month] (where month is longmonth not just months trapped), so can be used for phi~ndvi_1[i,m]. 





## array of dim [indiv,total months, max days] #  those not trapped=NA 
# matches up to p[i,m,d]
p.or.c <- p.or.c.array.fun(CH.secondary, covariate.data$temporal.covariates,n.sec.occ,list=TRUE)



### See model code in "robust_CJS_monthlylist_maxcov.bug"

################### Do all the data manipulation in R to create long data to remove loops in the observation code

obs.dat <- monthly.longdataCH.fun(CH.secondary, covariate.data$temporal.covariates, covariate.data$individual.covariates)


monthlyCH <- monthly.primaryCH.fun(CH.primary,covariate.data$temporal.covariates)


webmonths <- list()
for(i in 1:(length(session.list)-1)){
  x <- match(session.list[[i+1]],covariate.data$temporal.covariates$session)
  webmonths[[i]] <- x[which(is.finite(x))]
}
names(webmonths) <- names(session.list)[-1]

months.trapped.mat <- matrix(NA, nrow=dim(CH.secondary[[1]])[1],ncol=max(unlist(lapply(webmonths,length))))
length.months.trapped <- numeric()
for(i in 1:dim(months.trapped.mat)[1]){
  webnam <- covariate.data$individual.covariates$web[i] # this is a factor currently
  webi <- which(names(webmonths)==paste("web",webnam,sep="."))
  length.months.trapped[i]  <- length(webmonths[[webi]])
  months.trapped.mat[i,1:length.months.trapped[i]] <- webmonths[[webi]]
}

### need to make the code to different webs more general so can apply to other sites with different number and differently named webs
# same for model code (p loops)

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
  p.ind = replace(p.or.c,p.or.c==0,1),
  web = covariate.data$individual.covariates$web,
  sex = covariate.data$individual.covariates$sex,#
  n.sec.occ = n.sec.occ, # vector of number of secondary occasions 
  max.secondary.occasions = max(obs.dat$Sec), # not currently using
  n.primary.occasions = max(obs.dat$Prim), 
  monthlyCH = monthlyCH,
  z = known.state.cjs(monthlyCH),
  cjs.init.z=cjs.init.z,
  #CH.primary=CH.primary,
  n.months = dim(covariate.data$temporal.covariates)[1], 
  covariate.month = covariate.data$temporal.covariates$month, #1:12
  long.month = covariate.data$temporal.covariates$long.month, 
  months.trapped.mat = months.trapped.mat,
  length.months.trapped = length.months.trapped, 
  Prim = covariate.data$temporal.covariates$Prim, 
  ndvi_2 = covariate.data$ndvi_2,
  pt_0 = covariate.data$pt_0,
  pt_1 = covariate.data$pt_1,
  ShannonH_0 = covariate.data$ShannonH_0
  
  #lag.mat = matrix(c(covariate.data$ndvi_2,covariate.data$pt_0),ncol=2), # for the indicator variable
  #n.lags = 4 # for the indicator variable
) 


#initial values
inits=function(){list(z=cjs.init.z(monthlyCH,f.longmonth), mean.phi=runif(1,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1),alpha.0=runif(1,0,1),alpha.month=runif(11,0,1),  alpha.ndvi_2=runif(1,0,1), alpha.pt_0=runif(1,0,1), alpha.pt_1=runif(1,0,1), alpha.ShannonH_0=runif(1,0,1), alpha.male=runif(1,0,1),sigma.0=runif(1,0,1), sigma.recap=runif(1,0,1),sigma.male=runif(1,0,1),sigma.month=runif(11,0,1) )} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.month","alpha.ndvi_2", "alpha.pt_0","alpha.pt_1","alpha.ShannonH_0","alpha.male","sigma.0","sigma.recap","sigma.male","sigma.month")


date()
Z12.rCJS.diversity <- jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_monthlylist_diversity.bug",n.chains=3,n.thin=6,n.iter=10,n.burnin=5)
date() 
save.image("Z12rCJSdiversitymod.RData")

date()
Z12.rCJS.diversity2 <- jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_monthlylist_diversity2.bug",n.chains=2,n.thin=6,n.iter=10000,n.burnin=5000)
date() 
save.image("Z12rCJSdiversitymod.RData")




Z12.rCJS.maxcov
hist.plot.fun(Z12.rCJS.maxcov)
chain.plot.fun(Z12.rCJS.maxcov)
