###########################################################################
## Robust design Multistate: seronegative, seropositive, (unknown serostatus)
##  measuring survival (phi), recapture (p), probabilty of transition (Psi)
## For Zuni 1 and 2 starting at 199709. 
###########################################################################



#-----------------To Do------------------------------------------#
# 1) check known state and initial value functions
#     make sure the number categories match up
# 2) make sure obs.dat$State lines up with model code
# 3) Figure out how to deal wtih unknown serostatus
#  (Working with "RDMS_phi_sexinf.R" code)
#----------------------------------------------------------------#


## For some reason at Zuni2, few samples were tested for infection
##  until 199709. So I will start there.

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

load("UNMdata.RData")
sw.temp.data <- read.csv("sw_merge.csv")
source("CaptureHistoryFunctions.R")
source("JSMultistateInfFunctions.R")
source("RobustCJSfunctions.R")
source("DiversityFunctions.R")


session.list <- session.list.function(
  data=UNMcaptures, #assuming this data has all the dates (including those trapped but no animals/pema were caught) and is 'cleaned' and has a column called 'session' or 'Session' # pema.dirty in Emily's data
  site="Zuni",
  webs=c("1","2")
)

# remove sessions before 199709.
session.list <- lapply(session.list,function(x){x[which(x=="199709"):length(x)]})
  
MNAs.diversityZ12 <- diversity.df.function(
  data=UNMdata, # capture data 
  sites="Zuni", 
  webs=c("1","2"), 
  sessions=session.list$all.sessions, # sessions trapped 
  interpolate=FALSE, # do you want to interpolate NAs?
  scale=TRUE, # scale between 0 and 1 (divide by max)  
  include.pm = TRUE # reviewer wanted diversity indices calculated with pm removed (include.pm=FALSE), but removing pm leads to Inf in invSimpsonD calculations (when no other species are present)
)

UNMcaptures.cut <- UNMcaptures[order(UNMcaptures$Session),]
UNMcaptures.cut <- UNMcaptures.cut[-which(UNMcaptures.cut$Session<199709),]
CH.secondary <- MS.capture.history.function(
  data=UNMcaptures.cut, #assuming this data has all the dates (including those trapped but no animals/pema were caught) and is 'cleaned'
  site="Zuni",
  webs=c("1","2"),
  species="PM",
  SNV.unknown.state=TRUE # 1 = SNV neg,  2 = SNV positive, 3=unknown/not tested
)  
  

covariate.data <- monthly.covariate.function(
  capture.data=Zuni12.pema.data, #just pema so no overlap w other sites/species
  CH.secondary=CH.secondary, # list
  tags=rownames(CH.secondary[[1]]),
  sessions=session.list$all.sessions, #sessions trapped even if no pema caught
  temporal.data= sw.temp.data, # all the temporal covariate data in long format
  diversity.data= MNAs.diversityZ12,
  site= "Zuni", 
  web= c(1,2), 
  cov.list= list(pm=0,MNI=0)# this means I want MNI and MNA with no lag for each grid
)

CH.primary <- primary.MSch.fun(CH.secondary)

#rename secondary occasions from 1 to number of days 
for(i in 1:length(CH.secondary)){
  colnames(CH.secondary[[i]]) <- 1:dim(CH.secondary[[i]])[2]
}

# number of secondary occasions
n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))

monthlyCH <- monthly.primaryCH.fun(CH.primary,covariate.data$temporal.covariates)
# this works for multi-state

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


obs.dat <- monthly.longdataCH.fun(CH.secondary, covariate.data$temporal.covariates, covariate.data$individual.covariates)
# State = 0 for not caught
#       = 1 for caught as seronegative
#       = 2 for caught as seropositive
#       = 3 for caught as unknown serostatus
# This doesn't currently match up with the model code ("RDMS_phi_sexinf.R"). 
# Need to update one or both.


##### Bundle data
bugs.data <- list(
  nind = dplyr::n_distinct(obs.dat$ID), 
  n.obs = length(obs.dat$State),
  y = obs.dat$State, #observed as S or I or not observed
  prim = obs.dat$Prim, #primary occasion of observation 
  longmonth.obs = obs.dat$long.month, # longmonth of observation
  sec = obs.dat$Sec, #secondary occasion of obs
  id = obs.dat$ID, # individual (row of CH)
  f.longmonth = covariate.data$individual.covariates$f.longmonth, # longmonth first caught, for simulating z
  #p.or.c = p.or.c, # an array [i,m,d]to fit into p models
  #p.ind = replace(p.or.c,p.or.c==0,1),
  #web = covariate.data$individual.covariates$web,
  #sex = covariate.data$individual.covariates$sex,#
  I.dat = covariate.data$MNI_0, # needs to be a matrix : I.dat[i,t] like other temporal covariates
  N.dat = covariate.data$pm_0,
  n.sec.occ = n.sec.occ, # vector of number of secondary occasions 
  max.secondary.occasions = max(obs.dat$Sec), # not currently using
  n.primary.occasions = max(obs.dat$Prim), 
  monthlyCH = monthlyCH,
  z = known.state.SImsJS(monthlyCH), #change to ms
  jsmsinf.init=jsmsinf.init, 
  n.months = dim(covariate.data$temporal.covariates)[1], 
  long.month = covariate.data$temporal.covariates$long.month, 
  Prim = covariate.data$temporal.covariates$Prim 
) 

jsmsinf.init(ch=monthlyCH, num.aug=0)
#initial values
inits=function(){list(z=jsmsinf.init(ch=monthlyCH, num.aug=0), mean.p=runif(1,0,1),alpha.0=runif(1,0,1),alpha.t=runif(n.primary.occasions,0,1), alpha.maleS=runif(1,0,1),alpha.maleI=runif(1,0,1),alpha.inf=runif(1,0,1),sex.beta=runif(2,0,1) )} 

#parameters monitored
parameters=c("mean.p","alpha.0", "alpha.maleS","alpha.maleI","alpha.inf","sex.beta")




date()
Z12.rMS.phi.sexinf <- jags.parallel(data=bugs.data,inits,parameters,"RDMS_phi_sexinf.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() #
save.image("Z12rMSphisexinfmod.RData")

library(mcmcplots)
mcmcplot( Z12.rMS.phi.sexinf)

