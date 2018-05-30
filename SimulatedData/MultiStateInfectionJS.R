### Here, first I'm just making it so that can't go from Infected back to Susceptible
### (when put data in, need to make sure that's true in the data)

#############################################################
#
# Multistate Infection Joly-Seber capture-recapture model
#
##############################################################

library(R2jags)
setwd("~/Documents/JAGS")


# Estimation of survival (phi), force of infection (Psi), recruitment rates (f), S, I, & N population estimates


#############################################################
# 9.2.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiS <- 0.8             #survivial of Susceptible
phiI <- 0.7             #survival of Infected
psiSI <- 0.3            # prob of trasition from S to I
#psiIS <- 0              # prob of transisiton from I to S (should prob just get rid of this in process model)
pS <- 0.5               # capture prob of S
pI <- 0.6               # capture prob of I
n.occasions <- 6
n.states <- 3           #death is a state here, so S, I, & dead
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)  
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state at time t
# Dimension 2: state at time t+1
# Dimension 3: individual
# Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)      # number unique individuals [this is from the code file but the book says * (N.occasions-1)] but seems like it just sum(marked) is right
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1)) # a transition matrix for each individual (i) over time (t)
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiS*(1-psiSI), phiS*psiSI,     1-phiS,
      0,              phiI,           1-phiI,      #can't transition from I to S
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){              
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
#seen in S, I, not seen
        pS, 0,  1-pS,               # actually in S
        0,  pI, 1-pI,               # actually in I
        0,  0,  1                   # actually dead
        ), nrow = n.states, byrow = TRUE)
  } #t
} #i



# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}

# Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH


########################################################################## analyze the data

#  Analysis of the JS model as a multistate model
# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
rCH <- CH.ms          # Recoded CH
rCH[rCH==0] <- 3


#########################################################################
# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("msSIJS.bug")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phiS: survival probability of susceptible (S)
    # phiI: survival probability of infected (I)
    # psiSI:  probability of becoming infected (from S to I)
    #     can't go from I to S
    # pS: recapture probability as S
    # pI: recapture probability as I
    # gammaS: prob of entry into S class 
    # gammaI: prob of entry into I class (infected immigration)
    # -------------------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive as S
    # 3 alive as I
    # 4 dead
    # Observations (O):  
    # 1 seen as S 
    # 2 seen as I
    # 3 not seen
    # -------------------------------------------------
    
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      phiS[t] <- mean.phi[1]
      phiI[t] <- mean.phi[2]
      psiSI[t] <- mean.psi
      pS[t] <- mean.p[1]
      pI[t] <- mean.p[2]
      gammaS[t] ~ dunif(0, 1) #Prior for entry probabilities - why is this diff from p?
      gammaI[t] ~ dunif(0, 1) 
      #gammaS[t] <- mean.gamma[1]
      #gammaI[t] <- mean.gamma[2]
    }
    for (u in 1:2){     #for both states
      mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
      mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
      #mean.gamma[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
    }
    mean.psi ~ dunif(0, 1)    # Prior for mean transition

    # Define state-transition and observation matrices
    for (i in 1:M){  
      # Define probabilities of State(t+1) given State(t)
      for (t in 1:(n.occasions-1)){
      ### for proper process model, here say psiSI[t]<-beta*I[t] and define beta in priors instead of psiSI (need a logit or mlogit transform probably? Or is it fine given the priors?)
      ps[1,i,t,1] <- 1-gammaS[t]-gammaI[t]  # still not yet entered
      ps[1,i,t,2] <- gammaS[t]              # just entered as S
      ps[1,i,t,3] <- gammaI[t]              # just entered as I
      ps[1,i,t,4] <- 0                      # not yet entered to dead
      ps[2,i,t,1] <- 0                      # S to not yet entered
      ps[2,i,t,2] <- phiS[t] * (1-psiSI[t]) # survival of S to S
      ps[2,i,t,3] <- phiS[t] * psiSI[t]     # survival and transition from S to I
      ps[2,i,t,4] <- 1-phiS[t]              # S to dead
      ps[3,i,t,1] <- 0                      # I to not yet entered
      ps[3,i,t,2] <- 0                      # I to S (can't happen, tho check the data)
      ps[3,i,t,3] <- phiI[t]                # survival of I to I
      ps[3,i,t,4] <- 1-phiI[t]              # I to dead
      ps[4,i,t,1] <- 0                      # dead to not yet entered
      ps[4,i,t,2] <- 0                      # dead to S
      ps[4,i,t,3] <- 0                      # dead to I
      ps[4,i,t,4] <- 1                      # dead stay dead
    
      # Define probabilities of Observation(t) given State(t)
      # first index is state, last index is observation
      # could potentially include observation as wrong state (false neg or pos)
      po[1,i,t,1] <- 0        # not yet entered and observed as S
      po[1,i,t,2] <- 0        # not yet entered and observed as I
      po[1,i,t,3] <- 1        # not yet entered and not observed
      po[2,i,t,1] <- pS[t]    # in S and observed as S
      po[2,i,t,2] <- 0        # in S and observed as I 
      po[2,i,t,3] <- 1-pS[t]  # in S and not observed 
      po[3,i,t,1] <- 0        # in I and observed as S
      po[3,i,t,2] <- pI[t]    # in I and observed as I
      po[3,i,t,3] <- 1-pI[t]  # in I and not observed
      po[4,i,t,1] <- 0        # dead and observed as S
      po[4,i,t,2] <- 0        # dead and observed as I
      po[4,i,t,3] <- 1        # dead and not observed
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:M){
      # Define latent state at first occasion
      z[i,1] <- 1   # Make sure that all M individuals are in state 1 at t=1
      # No one has entered yet (state 1) at t=1, because when input data above (Ch.du), add a row before the actual data (where no has entered yet)
    
     for (t in 2:n.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1,]) #why t-1 here?
      } #t
    } #i

    # Calculate derived population parameters
    for (t in 1:(n.occasions-1)){
       qgamma[t] <- 1-gammaS[t]-gammaI[t]  # prob of not entering
    }
    cprob[1] <- gammaS[1]+gammaI[1]   # cummulative prob of entering (either as S or I)
    for (t in 2:(n.occasions-1)){
       cprob[t] <- (gammaS[t]+gammaI[t]) * prod(qgamma[1:(t-1)])
    } #t
    psi <- sum(cprob[])            # Inclusion probability
    for (t in 1:(n.occasions-1)){
       b[t] <- cprob[t] / psi      # Entry probability
    } #t
    
    for (i in 1:M){
       for (t in 2:n.occasions){
           Ss[i,t-1] <- equals(z[i,t], 2) #
           } #t
       for (t in 2:n.occasions){
           Is[i,t-1] <- equals(z[i,t], 3)
           } #t
       for (t in 1:(n.occasions-1)){
           dS[i,t] <- equals(z[i,t]-Ss[i,t],0) 
           dI[i,t] <- equals(z[i,t]-Is[i,t],0)
           } #t   
       aliveS[i] <- sum(Ss[i,])
       aliveI[i] <- sum(Is[i,])
     } #i
    
    for (t in 1:(n.occasions-1)){
       S[t]  <- sum(Ss[,t])       # Actual population size of S
       I[t]  <- sum(Is[,t])       # Actual pop size of I
       N[t]  <- S[t] + I[t]       # Actual total pop size
       BS[t] <- sum(dS[,t])       # Number of S entries
       BI[t] <- sum(dI[,t])       # Number of I entries
       B[t]  <- BS[t] + BI[t]     # total number of entries
       fS[t] <- BS[t]/N[t]        # per capita recruitment rate of S
       fI[t] <- BI[t]/N[t]        # per capita recruitment rate of I
       f[t]  <- B[t]/N[t]         # total per capita recruitment rate 
    } #t
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
bugs.data <- list(y = rCH, n.occasions = dim(rCH)[2], M = dim(rCH)[1], z = known.state.SImsJS(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(1, 0, 1), mean.p = runif(2, 0, 1), #mean.gamma = runif(2, 0, 1), 
                         z = SImsJS.init.z2(rCH))}  

# Parameters monitored
parameters <- c("mean.p", "mean.psi", "mean.phi", "Nsuper", "N", "f")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# run model 
date()
ms <- jags(bugs.data, inits, parameters, "msSIJS.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
date() 

#Error in jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
#Error in node z[201,3]
#Node inconsistent with parents


print(ms, digits = 3)






