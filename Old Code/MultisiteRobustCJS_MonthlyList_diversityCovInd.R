#################################################################
# test on smaller dataset

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
source("CaptureHistoryFunctions.R")

CH.secondary <- multisite.CJS.capture.history.function(dirty.data = pema.dirty, cleaned.data = pema.final.clean, site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"), species="PM")
  
session.list <- multisite.session.list.function(dirty.data = pema.dirty,
    site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2") )

#load("DiversityDataZuni12.RData") #or:
# source("DiversityFunctions.R")
# MNAs.diversityZ12 <- diversity.df.function(webs=1:2,interpolate=TRUE,scale=TRUE)

sw.temp.data <- read.csv("~/Documents/JAGS/BayesianMarkRecapSNV/EmilyThesisResearch/Data/norm_sw_merge.csv") 
CH.primary <- primary.ch.fun(CH.secondary)



# number of secondary occasions
n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))


covariate.data <- multisite.monthly.covariate.function(
  cleaned.data = pema.final.clean, # cleaned data
  CH.secondary = CH.secondary, # as monthly list
  tags = rownames(CH.secondary[[1]]), # tag names that line up to CH.secondary
  by.sitetag=TRUE, # if TRUE, then tags are specified by site e.g., "Zuni.1101"
  sessions = session.list$all.sessions, # all sessions to include e.g. 199806 (sessions trapped even if no pm caught)
  temporal.data= sw.temp.data, # data frame of monthly temporal data, with either a column 
  # called date or yearmon (must include months not trapped), 
  # must be in long format and have all the sites and all temporal data with no time lags
  site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"), 
  cov.list=list(ndvi=8,prcp=0,tmax=0) # list of temporal covariates and their time lags, 
  # this means use ndvi[t-8], prcp and tmax no lag 
)
# output is a list with the first element individual.covariates data frame
# second element called temporal.covariates to be used for long capture histories below
# other elements of the list are matrices of temporal covariates matching up individuals by web to their temporal data, e.g., temporal.covariates$ndvi_2 is a matrix of NDVI lag 2 with dimension [individual, month] (where month is longmonth not just months trapped), so can be used for phi~ndvi_1[i,m]. 


### stopped here


## array of dim [indiv,total months, max days] #  those not trapped=NA 
# matches up to p[i,m,d]
p.or.c <- multisite.p.or.c.array.fun(CH.secondary, covariate.data$temporal.covariates)



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


# put all the temporal covariates into an array for model selection, then will estimate an
# indicator which will give prob that covariate should be included in the model
n.covariates = 11
covariate.array = array(c(
  covariate.data$ndvi_2, #1
  covariate.data$pb_0, #2
  covariate.data$pf_0, #3
  covariate.data$pm_0, #4
  covariate.data$pt_0, #5
  covariate.data$pv_0, #6
  covariate.data$rg_0, #7
  covariate.data$ShannonH_0, #8
  covariate.data$speciesN_0, #9
  covariate.data$peros_0, #10
  covariate.data$other.sp_0 #11
  ),dim=c(dim(covariate.data$ndvi_2) , n.covariates))



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
  web = as.numeric(covariate.data$individual.covariates$web),
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
  covariate.array = covariate.array,
  n.covariates = n.covariates
) 


#initial values
inits=function(){list(z=cjs.init.z(monthlyCH,f.longmonth), mean.phi=runif(1,0,1),mean.p=runif(1,0,1),
  mean.c=runif(1,0,1),alpha.0=runif(1,0,1),alpha.month=runif(11,0,1),  
  alpha.ndvi_2=runif(1,0,1), alpha.pt_0=runif(1,0,1), alpha.pt_1=runif(1,0,1), 
  alpha.ShannonH_0=runif(1,0,1), alpha.male=runif(1,0,1),sigma.0=runif(1,0,1), 
  sigma.recap=runif(1,0,1),sigma.male=runif(1,0,1),sigma.month=runif(11,0,1),
  pind=runif(1,0,1), cov.coefT=runif(dim(covariate.array)[3],0,1) )} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.month","pind","cov.coefT","ind","alpha.male","sigma.0",
  "sigma.recap","sigma.male","sigma.month")


date()
Z12.rCJS.diversityCovInd <- jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_monthlylist_CovIndDiversity.bug",n.chains=3,n.thin=6,n.iter=15000,n.burnin=8000)
date() #36 hours
save.image("Z12rCJSdiversityCovInd.RData")




library(mcmcplots)
mcmcplot( Z12.rCJS.diversityCovInd)


Z12.rCJS.maxcov
hist.plot.fun(Z12.rCJS.maxcov)
chain.plot.fun(Z12.rCJS.maxcov)
