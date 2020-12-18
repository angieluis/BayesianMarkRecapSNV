#########################################################################################
#
# Estimation of survival, recruitment and population size using the Jolly-Seber model (multi-state formulation)
## Robust design with data as ragged array
#
##########################################################################################
library(R2jags)
library(dplyr)


# secondary capture history as a list of month matrices
CH.secondary <- 

# create primary capture histories
CH.primary <- primaryMS.fun(CH.secondary)

# Add dummy occasion
CH.primary.du <- primary.dummy.fun(CH.primary,notseen=3)

CH.secondary.du <- secondary.dummy.fun(CH.secondary,notseen=3)

# Augment data
num.aug <- 500
CH.primary.ms <- primary.augment.fun(CH.primary.du,notseen=3,num.aug=num.aug)
CH.secondary.ms <- secondary.augment.fun(CH.secondary.du,notseen=3,num.aug=num.aug)



source("JSRDMSAllConstant.R") #model code


################### Do all the data manipulation in R - create vectors
# and line up all the information and pass to BUGS with bugs.data

#  add the primary occassion to the data
for(i in 1:length(CH.secondary.ms)){
  CH.secondary.ms[[i]] <- cbind(data.frame(Prim = i), CH.secondary.ms[[i]])
}


obs_dat <- purrr::map_df(
  CH.secondary.ms,
  ~ tibble::as_tibble(.x) %>%    #
    dplyr::mutate(
      ID = 1:n()
    ) %>%
    tidyr::gather(Sec, Observation, -Prim, -ID) %>%
    dplyr::select(ID, Prim, Sec, Observation)
)


#  Subset observation data to observed bits (after first month)
y <- obs_dat[which(obs_dat$Prim>1),]


#### also need a column that indicates whether that individual
# has been caught before in that primary occasion (do we use p or c?)
#p.or.c <- numeric()
#for(i in 1:dim(y)[1]){
#  # the times that animal was caught that primary session
#  dat <- y[which(y$Prim==y$Prim[i] & y$ID==y$ID[i] & (y$Observation==1|y$Observation==2)),]

#  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
#    p.or.c[i] <- 0
#  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
#    firstcap <- min(as.numeric(dat$Sec))
#    p.or.c[i] <- ifelse(firstcap<y$Sec[i], 1 ,0)   # 0 is p, 1 is c
#  }
#}

#y <- data.frame(y,p.or.c)


# Bundle data
jags.data <- list(
  n.primary.occasions = max(y$Prim),
  # max.secondary.occasions = max(y$Sec),
  nind = dplyr::n_distinct(y$ID), # number of individuals (rows) including real and augmented
  y =  y$Observation, # Observation of animal (1 seen as S, 2 seen as I, 3 not seen)
  prim = y$Prim, # primary occasion
  # sec = y$Sec, # secondary occasion
  id = y$ID, # individual (1:nind)
  # p.or.c = y$p.or.c, # 0 if not caught before in that session, otherwise 1
  n.obs = nrow(y),
  z = known.state.SImsJS(ms=CH.primary.ms),
  jsmsinf.init = jsmsinf.init,
  CH.primary.ms = CH.primary.ms,
  num.aug = num.aug
)



inits <- function(){list(mean.phi = runif(2, 0, 1),
                         mean.p = runif(2, 0, 1),
                         mean.beta = runif(1, 0, 1),
                         # mean.c = runif(2, 0, 1),
                         # put in small values for gammas (because breaks o/w)
                         gammaS = rep(0.01, length(unique(y$Prim))),
                         gammaI = rep(0.01, length(unique(y$Prim))),
                         z = jsmsinf.init(CH.primary.ms, num.aug))}


# Parameters monitored
parameters <- c(
  "mean.p",
  # "mean.c",
  "mean.phi",
  "mean.beta",
  "b", "Nsuper", "N", "B", "f")



date()
js.rd.ms.inf <- jags.parallel(data=jags.data, inits, parameters, "js-rd-ms-inf.jags", n.chains = 3, n.thin = 6, n.iter = 10000, n.burnin = 5000)
date()

print(js.rd.ms.inf, digits = 3)

###
