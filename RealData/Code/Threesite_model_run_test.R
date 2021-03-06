########################################################################################
### Code to run combined data from Grand Canyon, Navajo, and Zuni together
########################################################################################


## This contains the covariates from the  Grand Canyon Max Model  ------------------- ##
 # where ndvi12, prcp6, tmin0, and tmax6 were selected as variables to include 


 # packages
   library(R2jags)
   library(tidyr)
   library(dplyr)


 # load RData specific to model run
   load("CombinedThreeSiteData.RData")

 # source functions and model
   source("Code/01_sw_data_functions.R")
   source("Code/combine_data_functions.R")
   source("Code/CombinedSites_GCModelCode.R")

   

   
## Bundle data -------------------------------------------------------------- ##


 # list of bugs data
   bugs.data <- list(
      nind                  = dplyr::n_distinct(obs.dat$ID),
      n.obs                 = length(obs.dat$State),
      y                     = obs.dat$State, 
      longmonth.obs         = obs.dat$long.month, 
      sec                   = as.numeric(obs.dat$Sec), # added as.numeric here
      id                    = obs.dat$ID,
      #ID.rm.last.day        = covariate.data$individual.covariates$ID[-which(covariate.data$individual.covariates$max.months==covariate.data$individual.covariates$f.longmonth)],  ## this is new, so don't model those caught on last day
      f.longmonth           = covariate.data$individual.covariates$f.longmonth, 
      p.or.c                = p.or.c,
      web                   = as.numeric(covariate.data$individual.covariates$web),
      site                  = as.numeric(factor(covariate.data$individual.covariates$site)),
      sex                   = as.numeric(covariate.data$individual.covariates$sex), 
      n.sec.occ             = n.sec.occ, 
      monthlyCH             = monthlyCH,
      z                     = known.state.cjs(monthlyCH),
      cjs.init.z.combined   = cjs.init.z.combined,
      n.months              = covariate.data$individual.covariates$max.months,  ## this varies by site so length=nind
      season                = covariate.data$season, # 1:4  ### now a matrix
      months.trapped.mat    = months.trapped.mat,
      length.months.trapped = length.months.trapped,
      Prim                  = covariate.data$Prim,  ## now a matrix
      ndvi                  = covariate.data$ndvi12_0,
      prcp                  = covariate.data$prcp6_0,
      tmin                  = covariate.data$tmin_0,
      tmax                  = covariate.data$tmax6_0,
      swe                   = covariate.data$swe_0,
      swe.winter            = covariate.data$swewinter_0
      
       )


## Initial values ----------------------------------------------------------- ##


 # supply initial values
   inits <- function() {
      list(
         z                  = bugs.data$cjs.init.z.combined(bugs.data$monthlyCH, bugs.data$f.longmonth), 
         alpha.0            = runif(1, 0, 1), 
         alpha.season       = runif(3, 0, 1), 
         alpha.web          = runif(max(bugs.data$web) - 1, 0, 1),
         alpha.male         = runif(1, 0, 1), 
         sigma.0            = runif(1, 0, 1),
         sigma.recap        = runif(1, 0, 1), 
         sigma.male         = runif(1, 0, 1), 
         sigma.season       = runif(3, 0, 1),
         sigma.web          = runif(max(bugs.data$web) - 1, 0, 1),
         alpha.swe          = runif(1, 0, 1),
         alpha.ndvi         = runif(1, 0, 1),
         alpha.prcp         = runif(1, 0, 1),
         alpha.tmin         = runif(1, 0, 1),
         alpha.tmax         = runif(1, 0, 1),
         alpha.ndvi.season  = runif(3, 0, 1),
         alpha.prcp.season  = runif(3, 0, 1),
         alpha.tmin.season  = runif(3, 0, 1),
         alpha.tmax.season  = runif(3, 0, 1),
         alpha.swe.winter   = runif(1, 0, 1)
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
                   "sigma.0", 
                   "sigma.recap", 
                   "sigma.male", 
                   "sigma.season"
                   
                    )


## Model specification ------------------------------------------------------ ##


 # date before run  
   date()


   Combined3sites.rCJS.GCModel <- jags(data     = bugs.data,
                                            inits, 
                                            parameters, 
                                            "Code/Combinedsites_GCModel.bug", 
                                            n.chains = 3, 
                                            n.thin   = 1, 
                                            n.iter   = 20, 
                                            n.burnin = 5)

   
   
   
   
   
   
   
   
   Combined3sites.rCJS.GCModel <- jags.parallel(data     = bugs.data,
                                                inits, 
                                                parameters, 
                                                "Code/Combinedsites_GCModel.bug", 
                                                n.chains = 3, 
                                                n.thin   = 6, 
                                                n.iter   = 30000, 
                                                n.burnin = 10000)
   
   
 # date after run
   date()


 # save session data
   save.image("Data/combined3sites_GCModelOutput.RData")

   


   
## -------------------------------------------------------------------------- ##
   
   
 # basic review
   
 
 # basic plots to check convergence and credible intervals,
   mcmcplot(Combined3sites.rCJS.GCModel,
            parms = c("alpha.0", 
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
                      "sigma.0", 
                      "sigma.recap", 
                      "sigma.male", 
                      "sigma.season"))
   # this takes awhile because of the derived parameters
   
   
 # MCMC vis package   
   MCMCplot(Combined3sites.rCJS.GCModel,
            xlim = c(-4, 4),
            ref_ovl = TRUE,
            params = c("alpha.0", 
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
                              "sigma.0", 
                              "sigma.recap", 
                              "sigma.male", 
                              "sigma.season")
            
              )
   