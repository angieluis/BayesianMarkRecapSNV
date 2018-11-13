#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## 
# 
# Data format: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)

#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")
#setwd("~/BayesianMarkRecapSNV/RealData")  # for PC

source("RobustCJSfunctions.R")
load("Z2pemaCH.RData")
#load("Z2temporalcovariates.RData")


CH.secondary <- Z2.pema.Ch.secondary
CH.primary <- primary.ch.fun(CH.secondary)
time.int <- Z2.primary.time.int.weeks 

### See model code in "robust_CJS_phi_month_p_dot_c_dot2.bug"

################### Do all the data manipulation in R - create vectors
# and line up all the information and pass to BUGS with bugs.data

#  A hack using a loop to add the primary occassion to the data
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
y <- dplyr::left_join(obs_dat, first_obs) %>%
  dplyr::filter(Prim >= f)


#### also need a column that indicates whether that individual
# has been caught before in that primary occasion (do we use p or c?)
p.or.c <- numeric()
for(i in 1:dim(y)[1]){
  # the times that animal was caught that primary session
  dat <- y[which(y$Prim==y$Prim[i] & y$ID==y$ID[i] & y$State==1),]
  
  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
    p.or.c[i] <- 0
  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
    firstcap <- min(as.numeric(dat$Sec))
    p.or.c[i] <- ifelse(firstcap<y$Sec[i], 1 ,0)   #BUGs doesn't like characters so 0 is p, 1 is c
  }
}

y <- data.frame(y,p.or.c)



##### Bundle data
bugs.data <- list(
  y = y$State,
  prim = y$Prim,
  sec = y$Sec,
  id = y$ID,
  f = first_obs$f, 
  p.or.c = y$p.or.c,
  nind = dplyr::n_distinct(y$ID), #n.secondary.occasions=n.secondary.occasions, 
  max.secondary.occasions = max(y$Sec), 
  n.primary.occasions = max(y$Prim), 
  time.int = time.int,
  n.obs = nrow(y),
  z = known.state.cjs(CH.primary),
  cjs.init.z=cjs.init.z,
  CH.primary=CH.primary,
  covariate.month = temporal.covariates$month,
  covariate.prim = temporal.covariates$covariate.prim,
  month.session = temporal.covariates$month.session
) 

#initial values
inits=function(){list(z=cjs.init.z(CH.primary,f),mean.phi=runif(12,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1))} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c")


date()
Z2.rcjs.phi.month.p.c.constant=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_phi_dot_p_dot_c_dot2.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() # 




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
