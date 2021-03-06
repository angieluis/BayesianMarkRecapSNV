#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
####### Data format: 
### Capture Histories: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)

###### Covariate dataframe called temporal.covariates
# temporal.covariates$long.month is every month of the start of the dataset 
  # to the end starting at 1 
# temporal.covariates$covariate.prim is each primary session. If there are
  # missing months in the dataset then there will be more than one long.month
  # to a covariate.prim (goes from 1 to the number of sessions sampled)
# temporal.covariates$month refers to month of the year (seasonality)

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
load("ZunipemaCH.RData")
#load("Zunitemporalcovariates.RData")
#load("Zuniindividualcovariates.RData")

CH.secondary <- Zuni.pema.Ch.secondary
CH.primary <- primary.ch.fun(CH.secondary)
time.int <- Zuni.primary.time.int.weeks 
temporal.covariates <- Zuni.temporal.covariates
individual.covariates <- Zuni.individual.covariates

### See model code in "robust_CJS_phi_month_p_dot_c_dot2.bug"

################### Do all the data manipulation in R - create vectors
# and line up all the information and pass to BUGS with bugs.data

#  add the primary occassion to the data
for(i in 1:length(CH.secondary)){
  CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
}


obs_dat <- purrr::map_df(
  CH.secondary, 
  ~ tibble::as_tibble(.x) %>%    # 
    dplyr::mutate(
      ID = 1:n()
    ) %>% 
    tidyr::gather(Sec, State, -Prim, -ID) %>%
    dplyr::select(ID, Prim, Sec, State)
)

first_obs <- obs_dat %>% 
  dplyr::group_by(ID) %>% 
  dplyr::summarise(f = min(Prim[State == 1]))

#  Subset observation data to observed bits
obs.dat <- dplyr::left_join(obs_dat, first_obs) %>%
  dplyr::filter(Prim >= f)

#### add individual covariates
obs.dat <- dplyr::left_join(obs.dat,individual.covariates)


#### also need a column that indicates whether that individual
# has been caught before in that primary occasion (do we use p or c?)
p.or.c <- numeric()
for(i in 1:dim(obs.dat)[1]){
  # the times that animal was caught that primary session
  dat <- obs.dat[which(obs.dat$Prim==obs.dat$Prim[i] & obs.dat$ID==obs.dat$ID[i] & obs.dat$State==1),]
  
  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
    p.or.c[i] <- 0
  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
    firstcap <- min(as.numeric(dat$Sec))
    p.or.c[i] <- ifelse(firstcap<obs.dat$Sec[i], 1 ,0)   #BUGs doesn't like characters so 0 is p, 1 is c
  }
}
p.or.c <- factor(p.or.c)



obs.dat <- dplyr::mutate(obs.dat, p.or.c)

#dummy data
#temporal.covariates$NDVI_1 <- rep(0.1,dim(temporal.covariates)[1])
#temporal.covariates$NDVI_2 <- rep(0.2,dim(temporal.covariates)[1])
#temporal.covariates$NDVI_3 <- rep(0.3,dim(temporal.covariates)[1])


#NDVI is the only temporal data that varies by individual/web
# so create a matrix of NDVI values by [i,m]
NDVI.m <- rbind(temporal.covariates$NDVI_1, temporal.covariates$NDVI_2, temporal.covariates$NDVI_3)
NDVI <- matrix(NA,ncol=dim(temporal.covariates)[1], nrow=dim(individual.covariates)[1])
for(i in 1:dim(individual.covariates)[1]){
  NDVI[i,] <- NDVI.m[individual.covariates$web[i],]
}


##### Bundle data
bugs.data <- list(
  y = obs.dat$State,
  prim = obs.dat$Prim,
  sec = obs.dat$Sec,
  id = obs.dat$ID,
  f = first_obs$f, 
  p.or.c = obs.dat$p.or.c,
  web = obs.dat$web,
  nind = dplyr::n_distinct(obs.dat$ID), 
  max.secondary.occasions = max(obs.dat$Sec), 
  n.primary.occasions = max(obs.dat$Prim), 
  time.int = time.int,
  n.obs = nrow(obs.dat),
  z = known.state.cjs(CH.primary),
  cjs.init.z=cjs.init.z,
  CH.primary=CH.primary,
  covariate.month = temporal.covariates$month,
  covariate.prim = temporal.covariates$covariate.prim,
  NDVI = NDVI,
  long.month = temporal.covariates$long.month
) 

#initial values
inits=function(){list(z=cjs.init.z(CH.primary,f),mean.phi=runif(12,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1),alpha.0=runif(1,0,1) ,alpha.1=runif(1,0,1))} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.1")


date()
Z2.rcjs.phi.month.p.c.constant=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_phi_month_p_dot_c_dot2.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() # 


date()
Z2.rcjs.phi.NDVI.p.c.constant=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_phi_NDVI_p_dot_c_dot.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() # 2 hours and 12 minutes


###############################################################################
### Model Comparisons
###############################################################################
mod.names <- objects()[grep("Z2.rcjs.",objects())]
mod.list <- list()
mod.table <- data.frame(model=mod.names,npar=rep(NA,length(mod.names)),
                        DIC=rep(NA,length(mod.names)))
for(i in 1:length(mod.names)){
  mod.list[[i]]<- get(mod.names[i]) 
  mod.table$DIC[i] <- mod.list[[i]]$BUGSoutput$DIC
  mod.table$npar[i] <- dim(mod.list[[i]]$BUGSoutput$summary)[1]-1
}
mod.table <- mod.table[order(mod.table$DIC),]
mod.table$delta.DIC <- mod.table$DIC-min(mod.table$DIC)
