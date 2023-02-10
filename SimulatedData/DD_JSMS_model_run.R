## Jolly Seber Multistate Model run code, Combined 3 sites using 
## Density Dependent K model - see "SimulateDensityDependence_noinf_array.R"
## Jolly-Seber models allow for individuals to enter the population
## from an augmented population - can estimate the prob of entrance & pop size


## To Do-------------------------------------------------------------------

## error with the ones trick. indexing?

## ------------------------------------------------------------------------

setwd("~/Documents/JAGS/BayesianMarkRecapSNV/SimulatedData")
load("SimDensDepArrayData.RData")
source("DD_JSMS_model_specification.R")
source("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/Code/01_sw_data_functions_more.R")

library(R2jags)

## for this version, the robust design capture histories are in an augmented array with dimensions
# CHarray.aug.du[i, m, d, w]
# CHarray[max number of individuals per web + number augmented, max number of 
# months per web + 1 for first dummy occasion, max number of secondary occasions per web, number of webs]
# for this simulation, 3 webs, all trapped all 40 months, and all had 3 secondary occasions
# the number of individuals varied by webs, so the number added to each web is

# number of individuals in 'original' capture data
n.inds
# number added to each web for data augmentation to get true pop size
n.add
#total individuals including real and augmented for each web (dimension of new array)
M <- n.inds + n.add 

# sex is an individual covariate and is a matrix of [i,w] with NAs for web 1 which didn't have the same n.inds
# need to give the 'fake' individuals a sex. Assume half are males and half are females
# first fill in the NAs, then add rows up to M
sex.aug <- sex
sind <- which(is.na(sex))
sex.aug[sind] <- rep(c(0,1),length=length(sind))
sex.aug <- rbind(sex.aug, matrix(c(0,1), nrow=min(n.add), ncol=dim(sex.aug)[2]))

## Need to figure out init function for this format

knownz <- array(NA, dim = dim(CHarray.aug.du)[-3])
for(w in 1:3){
  #convert array to list of monthly matrices so fits in function
  CHlist <- list()
  for(m in 1:dim(knownz)[2]){
    CHlist[[m]] <- CHarray.aug.du[ , m, , w]
  }
  monthlyCH <- primary.ch.fun(CHlist)  # as list of monthly matrices
  #function doesn't work for animals never observed (augmented)
  knownz[1:n.inds[w] ,  ,w] <- known.state.js(monthlyCH[1:n.inds[w], ]) 
  knownz[ ,1 ,w] <- NA # dummy occasion needs to be NA because in likelihood
  ## this so far 0=not yet entered and 1=alive. 
  ## But in the model state 1=not yet entered, 2=alive, 3=dead
  ## so add 1
  knownz[ , ,w] <- knownz[ , ,w] + 1
}

initz <- MSJSarray.init.z(ch=CHarray.aug.du, knownz)

jags.data <- list(
  n.webs = 3,
  n.inds =  M,# this should be length of webs and does it include augmented data?
  n.months = c(40, 40, 40)+1, # total months in z : length w +1 for dummy occasion
  length.months.trapped = c(40, 40, 40)+1, # how many months were trapped out of n.months : length w 
  months.trapped.mat = matrix(1:41, nrow=3, ncol=41, byrow=TRUE), # simplified here, but could be dimensions:[w, 1:length.months.trapped[w]]
  n.sec.occ = 3, # here simplified but could be dimensions:[w, Prim[i,m]]
  ndvi = ndvi, # covariate mat with dimensions: [w,m]
  season = season, # mat with dimensions: [w,m]
  sex = sex.aug, # mat with dimensions: [i,w]
  y = replace(CHarray.aug.du,CHarray.aug.du==0,2), # replace 0's with 2's 
  orig.y = CHarray.aug.du, # still with 0's for init function
  ones = matrix(1, nrow=3 ,ncol= dim(CHarray.aug.du)[2]), #[w,m] for the 'ones' trick
  knownz = knownz, # for the init function
  z = knownz, # for known states
  primary.ch.fun = primary.ch.fun, # for init function
  MSJSarray.init.z = MSJSarray.init.z
  
)

inits <- function(){list(m0 = runif(1, 0, 1), 
                         me = runif(1, 0, 1), 
                         m.male = runif(1, 0, 1), 
                         b0 = runif(1, 0, 1), # 
                         k.0 = runif(1, 0, 1), 
                         k.season = runif(3, 0, 1), # <- c(0.6,-0.4,-0.3) # spring, summer, fall, (winter is intercept)
                         k.ndvi = runif(1, 0, 1), # <- 0.35
                         sigma.0 = runif(1, 0, 1),
                         N1.est = rnorm(3, 20, 2), 
                         z = MSJSarray.init.z(orig.y, knownz)
                           )}


parameters <- c("m0",  #<- -2.944439 #logit(0.05)  #survival is 0.95 at N=0
              "me", # <- -0.8472979 #logit(0.3) # survival at equilibrium is 0.7
              "m.male", # <- 0.2 # additional mortality for males
              "b0", # <- 0.6931472 #log(2) # birth rate at N=0 is 2
              "be", # <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
              "k.0", # <- 3.1        # K is exp(k.0) at mean ndvi in winter 
              "k.season", # <- c(0.6,-0.4,-0.3) # spring, summer, fall, (winter is intercept)
              "k.ndvi", # <- 0.35
              "N1.est", # 
              "sigma.0", #            = 0.1 
              "N",
              "f")



date() #date and time before run
DD.JSMS.simulation.output <- jags(data     = jags.data,
                                            inits,
                                            parameters,
                                            "DD_JSMS_model_specification.bug",
                                            n.chains = 1,
                                            n.thin   = 0,
                                            n.iter   = 1,
                                            n.burnin = 0)

date() #date and time after run
save(DD.JSMS.simulation.output,file="DDJSMSsimoutput.RData")

# 'compilation phase' before initialization taking a long time. (didn't say 'initializing' or 'compiling' just 'glm module loaded')
# it is creating a graph representing the model in computer memory
# 1 day and 3.5 hours to run 1 iteration

#Graph information:
# Observed stochastic nodes: 252817
# Unobserved stochastic nodes: 81397
# Total graph size: 1078627

# Error in jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
#                       Error in node ones[1,2]
#                     Node inconsistent with parents