#########################################################################################
#
# Estimation of survival, recruitment and population size using the Jolly-Seber model (multi-state formulation)
# 
##########################################################################################
library(R2jags)
##################################
#  simulate data
##################################

# 10.4. Models with constant survival and time-dependent entry
# Define parameter values
n.occasions <- 7                         # Number of capture occasions
N <- 400                                 # Superpopulation size
phi <- rep(0.7, n.occasions-1)           # Survival probabilities
b <- c(0.34, rep(0.11, n.occasions-1))   # Entry probabilities 
p <- rep(0.5, n.occasions)               # Capture probabilities

PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions*N), ncol = n.occasions, nrow = N, byrow = T)

# Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI, P, b, N){
  B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
  n.occasions <- dim(PHI)[2] + 1
  CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
  # Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (t in 1:n.occasions){
    ent.occ <- c(ent.occ, rep(t, B[t]))
  }
  # Simulating survival
  for (i in 1:N){
    CH.sur[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
    if (ent.occ[i] == n.occasions) next
    for (t in (ent.occ[i]+1):n.occasions){
      # Bernoulli trial: has individual survived occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      ifelse (sur==1, CH.sur[i,t] <- 1, break)
    } #t
  } #i
  # Simulating capture
  for (i in 1:N){
    CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
  } #i
  # Full capture-recapture matrix
  CH <- CH.sur * CH.p
  
  # Remove individuals never captured
  cap.sum <- rowSums(CH)
  never <- which(cap.sum == 0)
  CH <- CH[-never,]
  Nt <- colSums(CH.sur)    # Actual population size
  return(list(CH=CH, B=B, N=Nt))
}

# Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH


##############################################
### Analyze data
##############################################

#  Analysis of the JS model as a multistate model
# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2                     # Not seen = 2, seen = 1


# 10.3.2. The JS model as a multistate model
# Specify model in BUGS language
sink("js-ms.jags")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # gamma: removal entry probability
    # p: capture probability
    #--------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    # Observations (O):
    # 1 seen 
    # 2 not seen
    #--------------------------------------
    
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      logit(phi[t]) <- mean.phi
      gamma[t] ~ dunif(0, 1) # Prior for entry probabilities - why is this diff from p? 
      logit(p[t]) <- mean.p
    }
    
    mean.phi ~ dnorm(0, 0.4)T(-10,10)    # Prior for mean survival
    mean.p ~ dnorm(0, 0.4)T(-10,10)      # Prior for mean capture
    

    # Define state-transition and observation matrices 	
    for (i in 1:M){  
    # Define probabilities of state S(t+1) given S(t)
      for (t in 1:(n.occasions-1)){
        ps[1,i,t,1] <- 1-gamma[t]    # not yet alive
        ps[1,i,t,2] <- gamma[t]     # entering (not yet alive to alive)
        ps[1,i,t,3] <- 0            # can't go from not yet alive to dead
        ps[2,i,t,1] <- 0            # can't go from alive to not yet alive
        ps[2,i,t,2] <- phi[t]       # prob of survival (alive to alive)
        ps[2,i,t,3] <- 1-phi[t]     # prob of dying (alive to dead)
        ps[3,i,t,1] <- 0            # can't go from dead to not yet alive
        ps[3,i,t,2] <- 0            # can't go from dead to alive
        ps[3,i,t,3] <- 1            # dead, stay dead
    
        # Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- 0            # prob of capture if not yet alive =0
        po[1,i,t,2] <- 1            # prob of not capture if not yet alive =1
        po[2,i,t,1] <- p[t]         # prob of capture if alive
        po[2,i,t,2] <- 1-p[t]       # prob of not capturing
        po[3,i,t,1] <- 0            # prob of capture if dead = 0 
        po[3,i,t,2] <- 1            # prob of not capturing if dead =1
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
        y[i,t] ~ dcat(po[z[i,t], i, t-1,]) #why t-1 here
      } #t
    } #i
    
    # Calculate derived population parameters
    for (t in 1:(n.occasions-1)){
      qgamma[t] <- 1-gamma[t]
    }
    cprob[1] <- gamma[1]
    for (t in 2:(n.occasions-1)){
      cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)]) # prob of anyone entering now, is prob of available (not yet entered) indiv entering now gamma[t] and prob of not entering at any time in the past, which is product of all previous (1-gamma)
    } #t

    psi <- sum(cprob[])            # Inclusion probability # the prob of entering the pop at any time step (even if never caught) - a proportion will never enter so not part of the dataset (one of the non-existent animals in the augmented dataset)

    for (t in 1:(n.occasions-1)){
      b[t] <- cprob[t] / psi      # Entry probability
    } #t
    
    for (i in 1:M){
      for (t in 2:n.occasions){
        al[i,t-1] <- equals(z[i,t], 2)             # why is this t-1? because starting with t=2 because t=1 was added? But don't do that for d below.
      } #t
      for (t in 1:(n.occasions-1)){
        d[i,t] <- equals(z[i,t]-al[i,t],0)
      } #t   
      alive[i] <- sum(al[i,])
    } #i
    
    for (t in 1:(n.occasions-1)){
      N[t] <- sum(al[,t])        # Actual population size
      B[t] <- sum(d[,t])         # Number of entries
      f[t] <- B[t]/N[t]         # recruitment rate  <----- I added this
    } #t
    for (i in 1:M){
      w[i] <- 1-equals(alive[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
    }
    ",fill = TRUE)
sink()


######################

# Bundle data
jags.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

# Initial values
# As always in JAGS, good initial values need to be specified for the latent state z. They need to correspond to the true state, which is not the same as the observed state. Thus, we have given initial values of "1" for the latent state at all places before an individual was observed, an initial value of "2" at all places when the individual was observed alive or known to be alive and an initial value of "3" at all places after the last observation. The folllowing function creates the initial values.
js.multistate.init <- function(ch, nz){
  ch[ch==2] <- NA
  state <- ch
  for (i in 1:nrow(ch)){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 2
  }
  state[state==0] <- NA
  get.first <- function(x) min(which(!is.na(x)))
  get.last <- function(x) max(which(!is.na(x)))   
  f <- apply(state, 1, get.first)
  l <- apply(state, 1, get.last)
  for (i in 1:nrow(ch)){
    state[i,1:f[i]] <- 1
    if(l[i]!=ncol(ch)) state[i, (l[i]+1):ncol(ch)] <- 3
    state[i, f[i]] <- 2
  }   
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  state[,1] <-NA ## I added this, and now it works (1st col is specified in likelihood)
  return(state)
} 

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = js.multistate.init(CH.du, nz))}    

# from the book and doesn't work:
# inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = cbind(rep(NA, dim(CH.ms)[1]), CH.ms[, -1]))}
  


# Parameters monitored
parameters <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B", "f")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3


date()
js.ms <- jags(jags.data, inits, parameters, "js-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
date() #working # took 1.5 hours, uh oh
# is it slow because no known-state function?


print(js.ms, digits = 3)

### Do I ignore the first time step in the output because it was added? OR is that already removed? (Seems to be already removed from some of the output)
