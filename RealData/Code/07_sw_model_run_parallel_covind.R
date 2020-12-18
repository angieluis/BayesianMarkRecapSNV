## Code to run multiple models on multiple cores ---------------------------- ##


 # load packages required
   library(R2jags)
   library(doParallel)
   library(foreach)
   library(rlist)

  
 # source function file - includes all packages needed to run models
   source("Code/01_sw_data_functions.R")

 
   
## Specify models sites and covariates -------------------------------------- ##
   
 
 # specify model names
   model.names <- c("Zuni.NDVI",
                    "Zuni.PRCP",
                    "Zuni.TMIN",
                    "Zuni.TMAX",
                    "Navajo.NDVI",
                    "Navajo.PRCP",
                    "Navajo.TMIN",
                    "Navajo.TMAX",
                    "GrandCanyon.NDVI",
                    "GrandCanyon.PRCP",
                    "GrandCanyon.TMIN",
                    "GrandCanyon.TMAX")
  
  
  # specify sites in same combination
    site <- c("Zuni",
              "Zuni",
              "Zuni",
              "Zuni",
              "Navajo",
              "Navajo",
              "Navajo",
              "Navajo",
              "GrandCanyon",
              "GrandCanyon",
              "GrandCanyon",
              "GrandCanyon")


 # specify covariates in same combination
   covs <- c("ndvi",
             "prcp",
             "tmin",
             "tmax",
             "ndvi",
             "prcp",
             "tmin",
             "tmax",
             "ndvi",
             "prcp",
             "tmin",
             "tmax") 
   
   
## Register the doParallel backend and specify foreach ---------------------- ##  
   
  
 # detect cores - 24
   detectCores() 
   
   
 # register cores - since we're running 12 models we'll register 12 cores
   registerDoParallel(cores = 12)
   
   
 # sanity check - how many workers is foreach going to use - 12
   getDoParWorkers()
   
   
 # time run
   timer <- proc.time()

      
## Specify foreach loops ---------------------------------------------------- ##
   
   
 # specify foreach that loops through:
 # loads capture histories, loads covariates
   models <- foreach(j = 1:length(model.names), .combine = c) %dopar% { 
     
     
     # NOTE - how to format or do we need nested foreach loops?
     # grab capture histories to feed next foreach loop 
     # load capture data (clean and dirty southwest data)
     load(paste("Data/Updated/", site[j], "_CaptureHistories.RData", sep = ""))
     
     
     # load covariate data (with capture histories)
     covariate.array <- get(paste(covs[j],".covariate.array",sep=""))
     
     
## List BUGS data ----------------------------------------------------------- ##
     
     
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
       
     season = covariate.data$temporal.covariates$season, # 1:4
     year   = as.numeric(factor(covariate.data$temporal.covariates$year)),
     months.trapped.mat    = months.trapped.mat,
     length.months.trapped = length.months.trapped,
     Prim                  = covariate.data$temporal.covariates$Prim,
     covariate.array       = covariate.array,
     n.covariates          = dim(covariate.array)[3]
     )
     
     
## Initial values ----------------------------------------------------------- ##
     
     
 # provides initial values (priors)
     inits <- function() {
       list(
         z = cjs.init.z(monthlyCH, bugs.data$f.longmonth), 
         alpha.0      = runif(1, 0, 1), 
         alpha.season = runif(3, 0, 1), 
         alpha.web    = runif(max(bugs.data$web) - 1, 0, 1),
         alpha.male   = runif(1, 0, 1), 
         sigma.0      = runif(1, 0, 1),
         sigma.recap  = runif(1, 0, 1),
         sigma.male   = runif(1, 0, 1), 
         sigma.season = runif(3, 0, 1),
         sigma.web    = runif(max(bugs.data$web) - 1, 0, 1),
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
   
     
 # run model  
 # generally...
 # n.chain = 3
 # n.thin  = 6
 # n.iter  = 15000 up to 20000 for second run
 # n.burn  = 8000 up to 12000 
   output <- jags(data = bugs.data,
                  inits, 
                  parameters,
                  "Code/CovInd_Interaction_Updated.bug",
                  n.chains = 3,
                  n.thin   = 6,
                  n.iter   = 50000,
                  n.burnin = 30000)
     
        
 # paste model names to output
   assign(paste(model.names[j],"output", sep = "_"), output)
     
     
 # save output
   save.image(paste("Data/ModelRuns/ModelRerun/Updated", 
                    model.names[j], 
                    "_CovInd_Updated.RData", sep = ""))
   names(output) <- paste(model.names[j], names(output), sep = ".")
   
   output
   
   }
   
   
 # get runtime
 # 50,000 iterations took 10.9 days
 # 20,000 iterations and 4 fewer variables took 3.8 days
   time.taken <- proc.time() - timer
   
   
 # save models
   save(models, file = "Data/ModelRuns/ModelRerun/Updated/models.RData")
   
## -------------------------------------------------------------------------- ##
   

 # NOTE: Poor convergence on alpha for this run (model run took approximately
 # five days at 20,000 iterations. Was unable to update the model run. Kept
 # receiving the error "model must be recompiled" because I'm a jackass and 
 # saved a JAGS model object in your workspace and then tried to use it in a new 
 # R session. Within the same R session the model object should remain valid and 
 # update code should have worked.
   
 # NOTE: Poor convergence on covariate models that have covariate inclusion.
 # This may be indicative of what Josh mentioned (see O'Hara and Silan.). 
 # 50,000 iterations took approximately 10 days.
   
   
## -------------------------------------------------------------------------- ##
   
   
 library(MCMCvis)
   
 # double-check
   zuni.ndvi <- models$Zuni.NDVI.BUGSoutput
   
   d
   
   
   