### Running a multistate infection model on simulated data


############################################################
## Problems:
##  still need to change model to add birth stuff
############################################################


# packages
library(R2jags)
library(tidyr)
library(dplyr)


# simulated data in "SimulateDensityDependence.R"
load("SimulatedData/SimDensDepMSinfData.RData")

# source functions and model
source("RealData/Code/01_sw_data_functions_more.R")
source("SimulatedData/DD_MSinf_model_specification.R")


## Bundle data -------------------------------------------------------------- ##


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
  I.dat                 = covariate.data$MNI_0,
  N.dat                 = covariate.data$pm_0
  
)


## Initial values ----------------------------------------------------------- ##


# supply initial values
inits <- function() {
  list(
    z                  = MSinf.init.z(monthlyCH,n.months), 
    m0                 = runif(1, 0, 1), 
    me                 = runif(1, 0, 1),
    m.male             = runif(1, 0, 1),
    m.infected         = runif(1, 0, 1),
    b0                 = runif(1, 0, 1),
    k.0                = runif(1, 0, 1),
    k.season           = runif(3, 0, 1), 
    k.ndvi             = runif(1, 0, 1),
    sigma.0            = runif(1, 0, 1),
    sigma.inf          = runif(1, 0, 1),
    rho                = runif(1, 0, 1),
    beta.male          = runif(1, 0, 1),
    beta.I             = runif(1, 0, 1)

    
  )
}


## Parameters monitored ----------------------------------------------------- ##


# parameters monitored
parameters <- c("m0", 
                "me",
                "m.male",
                "m.infected",
                "b0",
                "k.0",
                "k.season", 
                "k.ndvi",
                "sigma.0", 
                "sigma.inf",
                "rho",
                "beta.male",
                "beta.I"
)


# test.output <- jags(data     = bugs.data,
#                                                inits, 
#                                                parameters, 
#                                                "DD_MSinf_model_specification.bug", 
#                                                n.chains = 1, 
#                                                n.thin   = 2, 
#                                                n.iter   = 10, 
#                                                n.burnin = 4)

#Error in if (n.months[i] != max(n.months)) { : 
#    missing value where TRUE/FALSE needed
#  Called from: MSinf.init.z(monthlyCH, n.months)





#date before run
date()
DD.MSinf.simulation.output <- jags.parallel(data     = bugs.data,
                                    inits,
                                    parameters,
                                    "DD_MSinf_model_specification.bug",
                                    n.chains = 3,
                                    n.thin   = 5,
                                    n.iter   = 30000,
                                    n.burnin = 10000)

date() # 
save(DD.MSinf.simulation.output,file="SimulatedData/DDMSinfsimoutput.RData")


DD.MSinf.simulation.output

library(MCMCvis)
MCMCplot(DD.MSinf.simulation.output,xlim = c(-4, 4),ref_ovl = TRUE)




