### Here, first I'm just making it so that can't go from Infected back to Susceptible
### (when put data in, need to make sure that's true in the data)

#############################################################
#
# 9. Multistate capture-recapture models
#
##############################################################

library(R2jags)
setwd("~/Documents/JAGS")


# 9.2. Estimation of movement between two states


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

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


##################################################################
# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("msSI.bug")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phiS: survival probability of susceptible (S)
    # phiI: survival probability of infected (I)
    # psiSI:  probability of becoming infected (from S to I)
    # can't go from I to S
    # pS: recapture probability as S
    # pI: recapture probability as I
    # -------------------------------------------------
    # States (S):
    # 1 alive as S
    # 2 alive as I
    # 3 dead
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
    }
    for (u in 1:2){     #for both states
    mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
    }
    mean.psi ~ dunif(0, 1)    # Prior for mean transition

    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of State(t+1) given State(t)
    for (t in f[i]:(n.occasions-1)){
    ### for proper process model, here say psiSI[t]<-beta*I[t] and define beta in priors instead of psiSI (need a logit or mlogit transform probably? Or is it fine given the priors?)
    ps[1,i,t,1] <- phiS[t] * (1-psiSI[t])
    ps[1,i,t,2] <- phiS[t] * psiSI[t] 
    ps[1,i,t,3] <- 1-phiS[t]
    ps[2,i,t,1] <- 0              # because can't go from I to S (tho check the data)
    ps[2,i,t,2] <- phiI[t]        # ditto
    ps[2,i,t,3] <- 1-phiI[t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of Observation(t) given State(t)
    # could potentially include observation as wrong state (false neg or pos)
    po[1,i,t,1] <- pS[t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1-pS[t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pI[t]
    po[2,i,t,3] <- 1-pI[t]
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){

    # State process: draw State(t) given State(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])

    # Observation process: draw Obs(t) given State(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])

    } #t
    } #i
    }
    ",fill = TRUE)
sink()


# Function to create known latent states z 
###changed this so can't go back to S from I. 
#If observed as I, then not seen, and seen again later, when not seen must have been I. 
#If observed as S, not seen, then observed as S again, then must be S.
#Allows us to fill in a lot. And should speed up computation time
known.state.SIms=function(ms,notseen){ 
  # notseen: label for 'not seen' #here is 3 (Dead)
  state <- ms
  state[state==notseen] <- NA
  for(i in 1:dim(ms)[1]){
    if(length(which(ms[i,]==2))>0){ #filling in S's where can
    minI=min(which(ms[i,]==2))
    maxI=max(which(ms[i,]==2))
    state[i,minI:maxI]=2}
    if(length(which(ms[i,]==1))>0){  #filling in I's where can
      minI=min(which(ms[i,]==1))
      maxI=max(which(ms[i,]==1))
      state[i,minI:maxI]=1}
    state[i,min(which(ms[i,]<3))]=NA			
  }
  state[state==3]=NA
  return(state)
}





# Function to create initial values for unknown z  
#filling in something for every 3/NA (not captured) in known.state function after initial capture 
#randomly samples from alive states for those that are unknown after initial capture
# reversed from known.state - want NAs everywhere we know and 1:f. 
SIms.init.z <- function(ch, f){ 

  zch = known.state.SIms(ch, 3)
  states <- max(ch, na.rm = TRUE) #in current example, 3 states, 3rd is dead
  known.states <- 1:(states-1) # remove unknown state
  
  for (i in 1:dim(ch)[1]){ 
    v=which(is.na(zch[i,(f[i]+1):dim(zch)[2]])) # all unknown after first cap
    zch[i,-(v+f[i])] <- NA        # all the known states make NA #doesn't work when all time points are known because length of v is 0, so:
    rch=replace(ch,ch==3,NA)
    if(length(v)==0){zch[i,f[i]:dim(zch)[2]]<-NA}
    if(length(v)>0){
    for(j in 1:length(v)){
      zch[i,v[j]+f[i]] <- max(rch[i,1:(v[j]-1+f[i])],na.rm=TRUE) # for initial value gives max of those seen before (so won't return a 1 if after a 2)
      
    }      
    }
  }
  return(zch)
}






# Bundle data
bugs.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.SIms(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(1, 0, 1), mean.p = runif(2, 0, 1), z = SIms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# run model in jags  (~12 min)
date()
ms <- jags(bugs.data, inits, parameters, "msSI.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
date()



print(ms, digits = 3)






