## Model run ---------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##


 # MODEL CURRENTLY RUNNING
   # RERUN - ZUNI TMIN*SEASON 


 # MODELS FINISHED RUNNING
   # ZUNI NDVI*SEASON 8/28/2019
   # ZUNI PRCP*SEASON 8/30/2019
   # ZUNI TMIN*SEASON 9/09/2019
   # ZUNI TMAX*SEASON 9/10/2019
   # ZUNI SWE*SEASON  9/11/2019
   # NAVAJO NDVI*SEASON 9/14/2019 
   # NAVAJO PRCP*SEASON 9/15/2019 
   # NAVAJO SWE*SEASON 9/16/2019 
   # NAVAJO TMIN*SEASON 9/17/2019
   # NAVAJO TMAX*SEASON 9/18/2019
   # GRANDCANYON SWE*SEASON 9/19/2019
   # GRANDCANYON NDVI*SEASON 9/20/2019
   # GRANDCANYON PRCP*SEASON 9/21/2019
   # GRANDCANYON TMIN*SEASON 9/22/2019 (finished in <12 hours, weird AC, and 
     # convergence issues)
   # GRANDCANYON TMAX*SEASON 9/23/2019
   # RERUN - NAVAJO TMIN*SEASON 9/26/2019 (took 43 hours)


 # MODELS TO RUN
   # NAVAJO TMAX*SEASON 
   # GRANDCANYON NDVI*SEASON 
   # GRANDCANYON TMIN*SEASON 
   # GRANDCANYON TMAX*SEASON


 # IMPORTANT   
   # load .RData session pertinent to model run
   # specify covariate.array for model run
   # i.e., covariate.array = swe.covariate.array etc.

## Specify which model to run ----------------------------------------------- ##
## -------------------------------------------------------------------------- ##


 # clear environment

 # source functions
   source("Code/01_sw_data_functions.R")
   
   
 # load RData specific to model run
   load("Data/CaptureHistories/GrandCanyon_CaptureHistories.RData")

 
   
 # specify covariate array for this model run
   covariate.array = tmin.covariate.array
   
     
## Bundle data -------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
 
   
 # list of bugs data
   bugs.data <- list(
      nind  = dplyr::n_distinct(obs.dat$ID),
      n.obs = length(obs.dat$State),
      
      # observed yes/no
      y = obs.dat$State, 
      
      # longmonth of observation
      longmonth.obs = obs.dat$long.month, 
      
      # secondary occasion of obs
      sec = obs.dat$Sec, 
      
      # individual (row of CH)
      id = obs.dat$ID, 
      
      # longmonth first caught, for simulating
      # includes months not trapped
      f.longmonth = covariate.data$individual.covariates$f.longmonth, 
      
      # an array [i,m,d]to fit into p models
      p.or.c = p.or.c, 
      web    = as.numeric(covariate.data$individual.covariates$web),
      
      # under the individual covariate dataframe
      sex = covariate.data$individual.covariates$sex, 
      
      # vector of number of secondary occasions
      n.sec.occ = n.sec.occ, 
      monthlyCH = monthlyCH,
      z         = known.state.cjs(monthlyCH),
      
      # function used below
      cjs.init.z = cjs.init.z, 
      n.months   = dim(covariate.data$temporal.covariates)[1],
      
      # below commented out - if you use month instead of season
      # covariate.month = covariate.data$temporal.covariates$month, #1:12
      season = covariate.data$temporal.covariates$season, # 1:4
      year   = as.numeric(factor(covariate.data$temporal.covariates$year)),
      
      # below commented out - can't remember
      # long.month = covariate.data$temporal.covariates$long.month,
      months.trapped.mat    = months.trapped.mat,
      length.months.trapped = length.months.trapped,
      Prim                  = covariate.data$temporal.covariates$Prim,
      covariate.array       = covariate.array,
      n.covariates          = dim(covariate.array)[3]
      )
   
     
## Initial values ----------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##      
     
     
 # provides initial values (priors)
   inits <- function() {
     list(z = cjs.init.z(monthlyCH, f.longmonth),
          alpha.0      = runif(1, 0, 1), 
          alpha.season = runif(3, 0, 1), 
          alpha.web    = runif(max(web) - 1, 0, 1),
          alpha.male   = runif(1, 0, 1), 
          sigma.0      = runif(1, 0, 1),
          # sigma.recap  = runif(1, 0, 1), '
          sigma.recap  = 1.45,
           # specified this value based on first run output
          sigma.male   = runif(1, 0, 1), 
          sigma.season = runif(3, 0, 1),
          sigma.web    = runif(max(web) - 1, 0, 1),
          # sigma.year   = runif(max(year) - 1, 0, 1),
          # ran sigma.year (first run) see mcmc results for 
          # explanation of removal
          # pind = runif(1, 0, 1),
          cov.coefT = matrix(runif(dim(covariate.array)[3] * 4, 0, 1),
                             nrow = dim(covariate.array)[3],
                             ncol = 4)
          )
     }
     
     
## Parameters monitored ----------------------------------------------------- ##
     
     
 # parameters monitored - we can add things that aren't listed
 # won't monitor initial values or parameters that aren't specified in the code
   parameters <- c("alpha.0",
                   "alpha.male", 
                   "alpha.season", 
                   "alpha.web",
                   "cov.coefT", 
                   "ind", 
                   "sigma.0", 
                   "sigma.recap", 
                   "sigma.male", 
                   "sigma.season") 
     
     
## Model specification ------------------------------------------------------ ##
     
     
 # specify file 
 # model specific to covariate indicators, NOT specific to any one covariate
 # should not have to change model specification code to run differrent covs
   
   # this seems like redundant code
   #source("Code/05_sw_model_specification_covind.R")
     
     
 # date prior to model run
   date()
     
     
 # run model
 # generally...
 # n.chain = 3
 # n.thin  = 6
 # n.iter  = 15000 up to 20000 for second run, upped to 50000 to deal with
   # convergence issues
 # n.burn  = 8000 up to 12000 
   TEST.GC.rCJS.TMINCovInd <- jags.parallel(data = bugs.data, 
                                      inits, 
                                      parameters, 
                                      "Code/CovInd_Interaction.bug", 
                                      n.chains = 1, 
                                      n.thin   = 2, 
                                      n.iter   = 5, 
                                      n.burnin = 3)
     
 # date after model run
   date()
   
   
 # save output of model run - rename depends on model run
   save.image("Data/ModelRuns/REDO_Zuni_RDCJS_CovInd_TMIN_Inter.RData")
     
     
## -------------------------------------------------------------------------- ##