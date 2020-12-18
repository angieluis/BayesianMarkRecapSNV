### Here, first I'm just making it so that can't go from Infected back to Susceptible
### (when put data in, need to make sure that's true in the data)

#############################################################
#
# Robust Multistate Infection Joly-Seber capture-recapture model
#
##############################################################

library(R2jags)
setwd("~/Documents/JAGS")


# Estimation of survival (phi), force of infection (Psi), recruitment rates (f), 
# from Robust Design: p (capture), c (recapture),  S, I, & N population estimates


#############################################################
# Need to create function to simulate data


# Data format: Data frame where
# each row is an observation - each time an individual was caught, like the rodent data
# e.g.  individual, ID1, was caught in the third month on days 2 and 3 and was susceptible (state 1) and then again on month 5 first day and was infected (state 2):
# $individual  $primary   $secondary    $state
#   ID1           3           2           1
#   ID1           3           3           1
#   ID1           5           1           2
#   ID2           2           2           1 

observed.data # 


################################################################ analyze the data
#  Analysis of the Robust design JS model as a multistate model

ids <- unique(observed.data$individual)
nind <- length(ids)

n.primary.occasions <- length(unique(observed.data$primary))
max.secondary.occasions <- length(unique(observed.data$secondary))
n.secondary.occasions <- numeric()
for(m in 1:n.primary.occasions){
  n.secondary.occasions[m] <- length(unique(observed.data$secondary[which(observed.data$primary==m)]))
}


CH.primary <- matrix(0, ncol=n.primary.occasions ,nrow=nind)
for(i in 1:dim(observed.data)[1]){
  CH.primary[observed.data$individual[i], observed.data$primary[i]] <- observed.data$state[i]
}


# Add dummy occasion (new month at the beginning of the list)
CH.primary.du <- cbind(rep(0,dim(CH.primary)[1]),CH.primary)


# Augment data
nz <- 500
CH.secondary.ms <- lapply(CH.secondary.du, function(x){rbind(x,matrix(0, ncol = dim(x)[2], nrow = nz))})


# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
rCH.secondary <- CH.secondary.ms          # Recoded CH
rCH.secondary[rCH.secondary==0] <- 3

# create a recoded primary CH from the secondary capture history:
x <- lapply(rCH.secondary,rowSums)
v1 <- unlist(x)
rCH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(rCH.secondary[[1]])[1], ncol=length(rCH.secondary)) 

n.secondary.occasions <- unlist(lapply(rCH.secondary,function(x){dim(x)[2]}))
M = dim(rCH.primary)[1]

#########################################################################
#  Analysis of the model
# Specify model in BUGS language
sink("RobustmsSIJS.bug")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    #   phiS: survival probability of susceptible (S)
    #   phiI: survival probability of infected (I)
    #   psiSI:  probability of becoming infected (from S to I)
    #     can't go from I to S
    #   pS: capture probability as S
    #   cS: recapture (within the same primary session) probability as S
    #   pI: capture probability as I
    #   cI: recapture (within the same primary session) probability as I
    #   gammaS: prob of entry into S class 
    #   gammaI: prob of entry into I class (infected immigration)
    # -------------------------------------------------
    # States (S):
    #   1 not yet entered
    #   2 alive as S
    #   3 alive as I
    #   4 dead
    # Observations (O):  
    #   1 seen as S 
    #   2 seen as I
    #   3 not seen
    # -------------------------------------------------
    
    ############### Priors and constraints

    for(i in 1:nind){
      for(m in f[i]:n.primary.occasions){  ### for p need every time for phi need -1
        # phi, psi, gamma, have only 2 dimensions [indiv, and primary occasions]
        phiS[i,m] <- mean.phiS        # could specify covariates here
        phiI[i,m] <- mean.phiI        
        psiSI[i,m] <- mean.psiSI      
        gammaS[i,m] <- mean.gammaS    
        gammaI[i,m] <- mean.gammaI    
  
        # p and c have 3 dimensions [[primary]][indiv, secondary]		
        for(d in 1:n.secondary.occasions){
          pS[[m]][i, d] <- mean.pS  # could specify covariates here
          cS[[m]][i, d] <- mean.cS  
          pSI[[m]][i, d] <- mean.pI  
          cI[[m]][i, d] <- mean.cI  
        } #d for days
      } #m for months
    } #i for individual
    

    mean.phiS ~ dunif(0, 1)    # Prior for mean S survival
    mean.phiI ~ dunif(0, 1)    # Prior for mean I survival
    mean.pS ~ dunif(0, 1)      # Prior for mean S capture
    mean.pI ~ dunif(0, 1)      # Prior for mean I capture
    mean.cS ~ dunif(0, 1)      # Prior for mean S recapture
    mean.cI ~ dunif(0, 1)      # Prior for mean I recapture
    mean.psiSI ~ dunif(0, 1)   # Prior for mean S to I transition
    mean.gammaS ~ dunif(0, 1)  # Prior for mean prob of entry for S
    mean.gammaI ~ dunif(0, 1)  # Prior for mean prob of entry for I

    # Define state-transition and observation matrices
    for (i in 1:M){  
    # Define probabilities of State(t+1) given State(t)
      for (m in 1:(n.occasions-1)){
    ### for proper process model, here say psiSI[m]<-beta*I[m] and define beta in priors instead of psiSI (need a logit or mlogit transform probably? Or is it fine given the priors?)
        ps[1,i,m,1] <- 1-gammaS[i,m]-gammaI[i,m]  # still not yet entered
        ps[1,i,m,2] <- gammaS[i,m]                # just entered as S
        ps[1,i,m,3] <- gammaI[i,m]                # just entered as I
        ps[1,i,m,4] <- 0                          # not yet entered to dead
        ps[2,i,m,1] <- 0                          # S to not yet entered
        ps[2,i,m,2] <- phiS[i,m] * (1-psiSI[i,m]) # survival of S to S
        ps[2,i,m,3] <- phiS[i,m] * psiSI[i,m]     # survival and transition from S to I
        ps[2,i,m,4] <- 1-phiS[i,m]                # S to dead
        ps[3,i,m,1] <- 0                          # I to not yet entered
        ps[3,i,m,2] <- 0                          # I to S (can't happen, check data)
        ps[3,i,m,3] <- phiI[i,m]                  # survival of I to I
        ps[3,i,m,4] <- 1-phiI[i,m]                # I to dead
        ps[4,i,m,1] <- 0                          # dead to not yet entered
        ps[4,i,m,2] <- 0                          # dead to S
        ps[4,i,m,3] <- 0                          # dead to I
        ps[4,i,m,4] <- 1                          # dead stay dead
    
        # Define probabilities of Observation(t) given State(t)
        # first index is state, last index is observation
        # could potentially include observation as wrong state (false neg or pos)
        for (d in 1:n.secondary.occasions){
          # define prob of observation. If it's the first day of the primary session, then use p. Otherwise use p if not caught previously that session or c if was caught that session. 

### this means that can't have separate loops over months- need to all be one big loop where everything is defined for that month, because here, I'm refering to y which hasn't been calculated yet.. Will that work?

          p.effS <- ifelse(d==1, pS[[m]][i, d], ifelse(sum(y[[m]][i, 1:(d-1)])==0, pS[[m]][i, d], cS[[m]][i, d])) 
          p.effI <- ifelse(d==1, pI[[m]][i, d], ifelse(sum(y[[m]][i, 1:(d-1)])==0, pI[[m]][i, d], cI[[m]][i, d]))
          po[1,i,m,1] <- 0                # not yet entered and observed as S
          po[1,i,m,2] <- 0                # not yet entered and observed as I
          po[1,i,m,3] <- 1                # not yet entered and not observed
          po[2,i,m,1] <- p.effS           # in S and observed as S
          po[2,i,m,2] <- 0                # in S and observed as I 
          po[2,i,m,3] <- 1-p.effS         # in S and not observed 
          po[3,i,m,1] <- 0                # in I and observed as S
          po[3,i,m,2] <- p.effI           # in I and observed as I 
          po[3,i,m,3] <- 1-p.effI         # in I and not observed  
          po[4,i,m,1] <- 0                # dead and observed as S
          po[4,i,m,2] <- 0                # dead and observed as I
          po[4,i,m,3] <- 1                # dead and not observed
          } # d
        } # m
    } #i

########## need to put in same loops above
    # Likelihood  
    for (i in 1:M){
      # Define latent state at first occasion
      # dimensions [individual, primary session (month)]
      z[i,1] <- 1   # Make sure that all M individuals are in state 1 at t=1
      # No one has entered yet (state 1) at t=1, because when input data above (Ch.du), add a row before the actual data (where no has entered yet)
    
      for (m in 2:n.primary.occasions){
        # State process: draw S(m) given S(m-1)
        z[i,m] ~ dcat(ps[z[i,m-1], i, m-1, ])
        
        # Observation process: draw O(m) given S(m)
        # y has dimensions [[primary (month)]][individual, secondary (day)]
        y[[m]][i,d] ~ dcat(po[z[i,m], i, m-1, ]) # why m-1 here?
      } #m
    } #i
      
    # Calculate derived population parameters
    for (m in 1:(n.primary.occasions-1)){
      qgamma[m] <- 1-gammaS[m]-gammaI[m]  # prob of not entering
    }
    cprob[1] <- gammaS[1]+gammaI[1]   # cummulative prob of entering (either as S or I)
    for (m in 2:(n.primary.occasions-1)){
      cprob[m] <- (gammaS[m]+gammaI[t]) * prod(qgamma[1:(m-1)])
    } #m
    omega <- sum(cprob[])            # Inclusion probability
    for (m in 1:(n.primary.occasions-1)){
      b[m] <- cprob[m] / omega      # Entry probability
    } #t
    
    for (i in 1:M){
      for (m in 2:n.primary.occasions){
        Ss[i,m-1] <- equals(z[i,m], 2) #
        Is[i,m-1] <- equals(z[i,m], 3)
      } #m
      for (m in 1:(n.primary.occasions-1)){
        dS[i,m] <- equals(z[i,m]-Ss[i,m],0) 
        dI[i,m] <- equals(z[i,m]-Is[i,m],0)
      } #m   
      aliveS[i] <- sum(Ss[i,])
      aliveI[i] <- sum(Is[i,])
    } #i
    
    for (m in 1:(n.primary.occasions-1)){
      S[m]  <- sum(Ss[,m])       # Actual population size of S
      I[m]  <- sum(Is[,m])       # Actual pop size of I
      N[m]  <- S[m] + I[m]       # Actual total pop size
      BS[m] <- sum(dS[,m])       # Number of S entries
      BI[m] <- sum(dI[,m])       # Number of I entries
      B[m]  <- BS[m] + BI[m]     # total number of entries
      fS[m] <- BS[m]/N[m]        # per capita recruitment rate of S
      fI[m] <- BI[m]/N[m]        # per capita recruitment rate of I
      f[m]  <- B[m]/N[m]         # total per capita recruitment rate 
    } #m
    for (i in 1:M){
      w[i] <- 1-equals(aliveS[i],0)-equals(aliveI[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
  
    }
    ",fill = TRUE)
sink()


# Function to create known latent states z 
### fill in all known but unobserved states (can't go back to S from I)
#If observed as I, then not seen, and seen again later, when not seen must have been I. 
#If observed as S, not seen, then observed as S again, then must be S.
# Remember these are now states not observations, so coded differently. 
# 1 is not yet entered
# 2 is S
# 3 is I
# 4 is dead
# only able to fill in 2's and 3's
# Allows us to fill in a lot. And should speed up computation time
known.state.SImsJS <- function(ms, notseen){ # ms is multistate capture history
  # notseen: label for 'not seen' #here is 3 (Dead)
  state <- ms
  state[state==notseen] <- NA
  for(i in 1:dim(ms)[1]){
    if(length(which(ms[i, ] == 2)) > 0){ #filling in I's where can
      minI <- min(which(ms[i, ] == 2)) #I's are observation 2
      maxI <- max(which(ms[i, ] == 2))
      state[i, minI:maxI] <- 3}         # I's are state 3
    if(length(which(ms[i, ]==1)) > 0){  #filling in S's where can
      minS <- min(which(ms[i, ] == 1))  # S's are observation 1
      maxS <- max(which(ms[i, ] == 1))
      state[i, minS:maxS] <- 2}         # S's are state 2
  }
  return(state)
}




# Function to create initial values for unknown z  
#filling in something for every 3/NA (not captured) in known.state function  
#for those that are unknown, before it was known alive, make it not yet entered (1), once seen, make it alive and in same state as last seen (2 if last seen as S or 3 if last seen as I)
#  want NAs everywhere we know. 
SImsJS.init.z <- function(ch,states=4){ #states is number of states 
  kch <- known.state.SImsJS(ch, 3)   # known states
  zch <- matrix(NA,nrow=dim(kch)[1],ncol=dim(kch)[2])
  #  known.states <- 2:(states-1) # remove unobserved states (1 & 4)
  
  for (i in 1:dim(ch)[1]){ 
    v <- which(is.na(kch[i,])) # all unknown 
    #doesn't work when all time points are known because length of v is 0, so:
    if(length(v)==0){zch[i, ] <- NA}
    if(length(v)>0){
      for(j in 1:length(v)){
        zch[i, v[j]] <- ifelse(sum(kch[i,1:(v[j]-1)],na.rm=TRUE)==0,1,max(kch[i,1:(v[j]-1)],na.rm=TRUE)) # for initial value gives max of those seen before (so won't return a 1 if after a 2), if never seen before, sum==0, so put a 1 
        zch[ , 1] <- NA # in likelihood all animals are in 1 at t=1, so make NA
      }      
    }
  }
  return(zch)
}





# same as above but make dead after last time seen
SImsJS.init.z2 <- function(ch,states=4){ #states is number of states 
  kch <- known.state.SImsJS(ch, 3)   # known states
  zch <- matrix(NA,nrow=dim(kch)[1],ncol=dim(kch)[2])
  #  known.states <- 2:(states-1) # remove unobserved states (1 & 4)
  
  for (i in 1:dim(ch)[1]){ 
    v <- which(is.na(kch[i,])) # all unknown 
    #doesn't work when all time points are known because length of v is 0, so:
    if(length(v)==0){zch[i, ] <- NA}
    if(length(v)>0){
      for(j in 1:length(v)){
        zch[i, v[j]] <- ifelse(sum(kch[i,1:(v[j]-1)],na.rm=TRUE)==0,1,max(kch[i,1:(v[j]-1)],na.rm=TRUE)) # for initial value gives max of those seen before (so won't return a 1 if after a 2), if never seen before, sum==0, so put a 1 
      } #j     
      # if not seen at last occasion and ever seen, after last time seen, make dead (4)
      if(dim(ch)[2] %in% v & length(v)!=7){
        c <- 1:dim(ch)[2]
        zch[i, (max(c[-v])+1):max(c)] <- 4
      }
    } 
  } #i
  zch[ , 1] <- NA # in likelihood all animals are in 1 at t=1, so make NA
  return(zch)
}






# Bundle data
bugs.data <- list(y = rCH.secondary, n.primary.occasions = dim(rCH.secondary[[1]])[2] , n.secondary.ocassions = n.secondary.oaccsions, M = dim(rCH.secondary[[1]])[1], z = known.state.SImsJS(rCH.primary, 3))

# Initial values
inits <- function(){list(mean.phiS = runif(1, 0, 1), mean.phiI = runif(1, 0, 1), mean.psiSI = runif(1, 0, 1), mean.pS = runif(1, 0, 1), mean.pI = runif(1, 0, 1), mean.cS = runif(1, 0, 1), mean.cI = runif(1, 0, 1), mean.gammaS = runif(1, 0, 1), mean.gammaI = runif(1, 0, 1), 
          z = SImsJS.init.z2(rCH.primary))}  

# Parameters monitored
parameters <- c("mean.pS", "mean.pI", "mean.cS", "mean.cI", "mean.psiSI", "mean.phiS", "mean.phiI", "Nsuper", "N", "f")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# run model 
date()
ms <- jags(bugs.data, inits, parameters, "RobustmsSIJS.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
date() 

#Error in jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
#Error in node z[201,3]
#Node inconsistent with parents


print(ms, digits = 3)






