### Running a multistate infection model on simulated data


############################################################
## Problems:
##  
############################################################


# packages
library(R2jags)
library(tidyr)
library(dplyr)


# simulated data in "MSinf_3site_SimulateData.R"
load("SimulatedData/CombinedSimMSDatasimpler2.RData")

# source functions and model
source("RealData/Code/01_sw_data_functions_more.R")
source("SimulatedData/MSinf_3site_model_specification.R")


## Bundle data -------------------------------------------------------------- ##

maxI <- 20

# list of bugs data
bugs.data <- list(
  nind                  = dplyr::n_distinct(obs.dat$ID),
  n.obs                 = length(obs.dat$State),
  y                     = obs.dat$State, # 1=observed as S; 2=observed as I; 3=not observed
  longmonth.obs         = obs.dat$long.month, 
  sec                   = as.numeric(obs.dat$Sec), # added as.numeric here
  id                    = obs.dat$ID,
  #ID.rm.last.day        = covariate.data$individual.covariates$ID[-which(covariate.data$individual.covariates$max.months==covariate.data$individual.covariates$f.longmonth)],  ## this is new, so don't model those caught on last day
  f.longmonth           = covariate.data$individual.covariates$f.longmonth, 
  f.state               = covariate.data$individual.covariates$f.state,
  p.or.c                = p.or.c,
  web                   = as.numeric(covariate.data$individual.covariates$web),
  site                  = as.numeric(factor(covariate.data$individual.covariates$site)),
  sex                   = as.numeric(covariate.data$individual.covariates$sex), 
  n.sec.occ             = n.sec.occ, 
  monthlyCH             = monthlyCH,
  z                     = known.state.SImsInf(monthlyCH),
  MSinf.init.z          = MSinf.init.z,
  known.state.SImsInf   = known.state.SImsInf,
  n.months              = covariate.data$individual.covariates$max.months,  ## this varies by site so length=nind
  season                = covariate.data$season, # 1:4  ### now a matrix
  months.trapped.mat    = months.trapped.mat,
  length.months.trapped = length.months.trapped,
  Prim                  = covariate.data$Prim,  ## now a matrix
  ndvi                  = covariate.data$ndvi_0,
  prcp                  = covariate.data$prcp_0,
  tmin                  = covariate.data$tmin_0,
  tmax                  = covariate.data$tmax_0,
  swe                   = covariate.data$swe_0,
  swe.winter            = covariate.data$swewinter_0,
  maxI                  = 20,
  I.dat                 = covariate.data$MNI_0/maxI
  
)


## Initial values ----------------------------------------------------------- ##


# supply initial values
inits <- function() {
  list(
    z                  = MSinf.init.z(monthlyCH,n.months), 
    alpha.0            = runif(1, 0, 1), 
    alpha.season       = runif(3, 0, 1), 
    alpha.ndvi         = runif(1, 0, 1),
    alpha.inf          = runif(1, 0, 1),
    alpha.inf.male     = runif(1, 0, 1),
    sigma.0            = runif(1, 0, 1),
    sigma.inf          = runif(1, 0, 1),
    beta.0             = runif(1, 0, 1),
    beta.male          = runif(1, 0, 1),
    beta.I             = runif(1, 0, 1)

    
  )
}


## Parameters monitored ----------------------------------------------------- ##


# parameters monitored
parameters <- c("alpha.0", 
                "alpha.season", 
                "alpha.ndvi",
                "alpha.inf ",
                "alpha.inf.male" ,
                "sigma.0", 
                "sigma.inf",
                "beta.0",
                "beta.male",
                "beta.I"
)


#test.output <- jags(data     = bugs.data,
#                                               inits, 
#                                               parameters, 
#                                               "MSinf_3site_model_specification.bug", 
#                                               n.chains = 1, 
#                                               n.thin   = 2, 
#                                               n.iter   = 10, 
#                                               n.burnin = 4)




# date before run  
date()
MSinf.3site.simulation.output <- jags.parallel(data     = bugs.data,
                                     inits, 
                                     parameters, 
                                     "MSinf_3site_model_specification.bug", 
                                     n.chains = 3, 
                                     n.thin   = 5, 
                                     n.iter   = 30000, 
                                     n.burnin = 10000)

date() # about 16 hours
save(MSinf.3site.simulation.output,file="SimulatedData/MSinf3sitesimoutput.RData")


MSinf.3site.simulation.output
#Inference for Bugs model at "MSinf_3site_model_specification.bug", fit using jags,
#3 chains, each with 30000 iterations (first 10000 discarded), n.thin = 5
#n.sims = 12000 iterations saved
#                  mu.vect sd.vect     2.5%      25%       50%       75%
#alpha.0             0.954   0.101    0.757    0.887     0.952     1.023
#alpha.inf          -0.951   0.155   -1.255   -1.055    -0.951    -0.847
#alpha.inf.male     -0.369   0.196   -0.751   -0.503    -0.372    -0.236
#alpha.ndvi          0.250   0.078    0.097    0.197     0.250     0.302
#alpha.season[1]     0.901   0.163    0.586    0.793     0.899     1.009
#alpha.season[2]    -0.511   0.150   -0.803   -0.613    -0.510    -0.410
#alpha.season[3]    -0.272   0.170   -0.603   -0.386    -0.273    -0.159
#beta.0             -1.853   0.155   -2.154   -1.958    -1.853    -1.747
#beta.I              1.424   0.445    0.544    1.128     1.426     1.721
#beta.male           0.497   0.154    0.196    0.392     0.495     0.601
#sigma.0             0.227   0.031    0.165    0.205     0.227     0.248
#sigma.inf           0.053   0.068   -0.077    0.006     0.053     0.099
#deviance        10015.811  47.689 9925.424 9983.377 10015.113 10047.853

                     #97.5%  Rhat n.eff
#alpha.0             1.155 1.001  7700
#alpha.inf          -0.652 1.001  4200
#alpha.inf.male      0.011 1.001  4800
#alpha.ndvi          0.403 1.001 12000
#alpha.season[1]     1.227 1.001  4400
#alpha.season[2]    -0.220 1.001 12000
#alpha.season[3]     0.060 1.001  8700
#beta.0             -1.551 1.001 12000
#beta.I              2.297 1.001 12000
#beta.male           0.804 1.001 12000
#sigma.0             0.289 1.001  8300
#sigma.inf           0.186 1.001  6900
#deviance        10109.760 1.001 12000

#For each parameter, n.eff is a crude measure of effective sample size,
#and Rhat is the potential scale reduction factor (at convergence, Rhat=1#).

#DIC info (using the rule, pD = var(deviance)/2)
#pD = 1137.2 and DIC = 11153.1
#DIC is an estimate of expected predictive error (lower deviance is better).

########################## not bad
# actual values:
#alpha.0            = 1 
#alpha.season       = c(1,-0.5,0) 
#alpha.ndvi         = 0.2
#alpha.inf          = -1
#alpha.inf.male     = -0.5
#sigma.0            = 0.1 
#sigma.inf          = 0.1
#beta.0             = -2
#beta.male          = 0.5
#beta.I             = 2

library(MCMCvis)
MCMCplot(MSinf.3site.simulation.output,xlim = c(-4, 4),ref_ovl = TRUE)




