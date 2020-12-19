


# packages
library(R2jags)
library(tidyr)
library(dplyr)


# simulated data in "SimulatedData/MultistateInf3steCovariates.R"
load("SimulatedData/CombinedSimMSData.RData")

# source functions and model
source("RealData/Code/01_sw_data_functions_more.R")
source("RealData/Code/Multistate_Combined_GCModelCode2.R")


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
  I.dat                 = covariate.data$MNI_0
  
)


## Initial values ----------------------------------------------------------- ##


# supply initial values
inits <- function() {
  list(
    z                  = bugs.data$MSinf.init.z(bugs.data$monthlyCH), 
    alpha.0            = runif(1, 0, 1), 
    alpha.season       = runif(3, 0, 1), 
    alpha.web          = runif(max(bugs.data$web) - 1, 0, 1),
    alpha.male         = runif(1, 0, 1), 
    alpha.swe          = runif(1, 0, 1),
    alpha.ndvi         = runif(1, 0, 1),
    alpha.prcp         = runif(1, 0, 1),
    alpha.tmin         = runif(1, 0, 1),
    alpha.tmax         = runif(1, 0, 1),
    alpha.ndvi.season  = runif(3, 0, 1),
    alpha.prcp.season  = runif(3, 0, 1),
    alpha.tmin.season  = runif(3, 0, 1),
    alpha.tmax.season  = runif(3, 0, 1),
    alpha.swe.winter   = runif(1, 0, 1),
    alpha.inf          = runif(1, 0, 1),
    alpha.inf.male     = runif(1, 0, 1),
    sigma.0            = runif(1, 0, 1),
    sigma.recap        = runif(1, 0, 1), 
    sigma.male         = runif(1, 0, 1), 
    sigma.season       = runif(3, 0, 1),
    sigma.web          = runif(max(bugs.data$web) - 1, 0, 1),
    sigma.inf          = runif(1, 0, 1),
    sigma.inf.male     = runif(1, 0, 1),
    beta.0             = runif(1, 0, 1),
    beta.male          = runif(1, 0, 1),
    beta.site          = runif(max(bugs.data$site) - 1, 0, 1),
    immig              = runif(1, 0, 1)
    
    
  )
}


## Parameters monitored ----------------------------------------------------- ##


# parameters monitored - we can add things that aren't listed
# won't monitor initial values or parameters that aren't specified in the code
# see Zuni model run for discussion on monitoring phi
parameters <- c("alpha.0", 
                "alpha.male", 
                "alpha.season", 
                "alpha.web", 
                "alpha.swe",
                "alpha.ndvi",
                "alpha.prcp",
                "alpha.tmin", 
                "alpha.tmax",
                "alpha.ndvi.season",
                "alpha.prcp.season",
                "alpha.tmin.season",
                "alpha.tmax.season",
                "alpha.swe.winter", 
                "alpha.inf ",
                "alpha.inf.male" ,
                
                "sigma.0", 
                "sigma.recap", 
                "sigma.male", 
                "sigma.season",
                "sigma.web",
                "sigma.inf",
                "sigma.inf.male",
                "beta.0",
                "beta.male" ,
                "beta.site" ,
                "immig"              
                
                
)


## Model specification ------------------------------------------------------ ##


# date before run  
date()


Combined3sites.MSinf.GCModel <- jags(data     = bugs.data,
                                     inits, 
                                     parameters, 
                                     "Multistate_Combinedsites_GCModel.bug", 
                                     n.chains = 3, 
                                     n.thin   = 1, 
                                     n.iter   = 20, 
                                     n.burnin = 5)





