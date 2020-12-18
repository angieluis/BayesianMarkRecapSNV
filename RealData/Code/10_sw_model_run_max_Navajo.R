## Navajo Max Model  -------------------------------------------------------- ##


 # code to run Navajo 1 and 2 robust design CJS model
 # where ndvi0, tmin0, and tmax0 were selected as variables to include 


 # packages
   library(R2jags)
   library(tidyr)
   library(dplyr)


## Specify which model to run ----------------------------------------------- ##
## -------------------------------------------------------------------------- ##


 # clear environment


 # source functions
   source("Code/01_sw_data_functions.R")


 # load RData specific to model run
   load("Data/Updated/Navajo_CaptureHistories_Jan2020.RData")


## Bundle data -------------------------------------------------------------- ##


 # list of bugs data
   bugs.data <- list(
      nind                  = dplyr::n_distinct(obs.dat$ID),
      n.obs                 = length(obs.dat$State),
      y                     = obs.dat$State, 
      longmonth.obs         = obs.dat$long.month, 
      sec                   = obs.dat$Sec, 
      id                    = obs.dat$ID,
      f.longmonth           = covariate.data$individual.covariates$f.longmonth, 
      p.or.c                = p.or.c,
      web                   = as.numeric(covariate.data$individual.covariates$web),
      sex                   = covariate.data$individual.covariates$sex, 
      n.sec.occ             = n.sec.occ, 
      monthlyCH             = monthlyCH,
      z                     = known.state.cjs(monthlyCH),
      cjs.init.z            = cjs.init.z,
      n.months              = dim(covariate.data$temporal.covariates)[1],
      season                = covariate.data$temporal.covariates$season, # 1:4
      months.trapped.mat    = months.trapped.mat,
      length.months.trapped = length.months.trapped,
      Prim                  = covariate.data$temporal.covariates$Prim,
      ndvi                  = covariate.data$ndvi_0,
      tmin                  = covariate.data$tmin_0,
      tmax                  = covariate.data$tmax_0,
      swe                   = covariate.data$swe_0,
      swe.winter            = covariate.data$swewinter_0,
      
      # derived parameters
      male.web1.index       = 149,
      male.web2.index       = 116,
      female.web1.index     = 117,
      female.web2.index     = 143
      
      )



## Initial values ----------------------------------------------------------- ##


 # supply initial values
   inits <- function() {
      list(
         z                  = cjs.init.z(monthlyCH, f.longmonth), 
         alpha.0            = runif(1, 0, 1), 
         alpha.season       = runif(3, 0, 1), 
         alpha.web          = runif(max(web) - 1, 0, 1),
         alpha.male         = runif(1, 0, 1), 
         sigma.0            = runif(1, 0, 1),
         sigma.recap        = runif(1, 0, 1), 
         sigma.male         = runif(1, 0, 1), 
         sigma.season       = runif(3, 0, 1),
         sigma.web          = runif(max(web) - 1, 0, 1),
         alpha.swe          = runif(1, 0, 1),
         # alpha.ndvi         = runif(1, 0, 1),
         # alpha.tmin         = runif(1, 0, 1),
         # alpha.tmax         = runif(1, 0, 1),
         alpha.ndvi.season  = runif(3, 0, 1),
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
                   # "alpha.ndvi",
                   # "alpha.tmin", 
                   # "alpha.tmax",
                   "alpha.ndvi.season",
                   "alpha.tmin.season",
                   "alpha.tmax.season",
                   "alpha.swe.winter", 
                   "sigma.0", 
                   "sigma.recap", 
                   "sigma.male", 
                   "sigma.season",
                   
                   # derived parameters
                   "phi.male.web1",
                   "phi.male.web2",   
                   "phi.female.web1", 
                   "phi.female.web2",
                   "p.male.web1",
                   "p.male.web2",
                   "p.female.web1",
                   "p.female.web2"
                   )


## Model specification ------------------------------------------------------ ##


 # date before run  
   date()


   N12.rCJS.MaxModel.J22 <- jags.parallel(data     = bugs.data,
                                          inits, 
                                          parameters, 
                                          "Code/Navajo_MaxModel_June2020.bug", 
                                          n.chains = 3, 
                                          n.thin   = 6, 
                                          n.iter   = 100000, 
                                          n.burnin = 40000)


 # date after run
   date()


 # save session data
   save.image("Data/ModelRuns/ForAnalysis/Navajo_MaxModel_June2020.RData")


## Update ------------------------------------------------------------------- ##
   
   
 # recompile - this has to happen each time a session ends
   recompile(N12.rCJS.MaxModel.F24,
             progress.bar = "text")
   
   
 # date before run
   date()
   
   
 # update - worked with default iterations and thinning (1000 and 1)
   update.N12.rCJS.MaxModel <- autojags(N12.rCJS.MaxModel.F24,
                                        n.iter = 50000,
                                        n.thin = 6,
                                        progress.bar = "text")
   
 # date after run
   date()
   
   
 # save session data
   save.image("Data/ModelRuns/ForAnalysis/Navajo_MaxModel_Mar2020.RData")
   
   
## -------------------------------------------------------------------------- ##
   
   
 # basic review
   
   
 # basic plots to check convergence and credible intervals
   mcmcplot(N12.rCJS.MaxModel.F24)
   # this takes awhile because of the derived parameters
   
   
 # MCMC vis package   
   MCMCplot(N12.rCJS.MaxModel.F24,
            xlim = c(-4, 4),
            ref_ovl = TRUE,
            params = c("alpha.0", 
                       "alpha.male", 
                       "alpha.season", 
                       "alpha.web", 
                       "alpha.swe", 
                       "alpha.tmin", 
                       "alpha.tmax",
                       "alpha.tmin.season",
                       "alpha.tmax.season",
                       "alpha.swe.winter", 
                       "sigma.0", 
                       "sigma.recap", 
                       "sigma.male", 
                       "sigma.season"))
   