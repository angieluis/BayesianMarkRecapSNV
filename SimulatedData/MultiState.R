#############################################################
#
# 9. Multistate capture-recapture models
#
##############################################################

# 9.2. Estimation of movement between two sites
# 9.2.1. Model description
# 9.2.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.8
phiB <- 0.7
psiAB <- 0.3
psiBA <- 0.5
pA <- 0.7
pB <- 0.4
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)  
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB), phiA*psiAB,     1-phiA,
      phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
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


#######################################################


# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("ms.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phiA: survival probability at site A
    # phiB: survival probability at site B
    # psiAB: movement probability from site A to site B
    # psiBA: movement probability from site B to site A
    # pA: recapture probability at site A
    # pB: recapture probability at site B
    # -------------------------------------------------
    # States (S):
    # 1 alive at A
    # 2 alive at B
    # 3 dead
    # Observations (O):  
    # 1 seen at A 
    # 2 seen at B
    # 3 not seen
    # -------------------------------------------------
    
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
    }
    for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
    mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
    }
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
    ps[1,i,t,2] <- phiA[t] * psiAB[t]
    ps[1,i,t,3] <- 1-phiA[t]
    ps[2,i,t,1] <- phiB[t] * psiBA[t]
    ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
    ps[2,i,t,3] <- 1-phiB[t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- pA[t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1-pA[t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pB[t]
    po[2,i,t,3] <- 1-pB[t]
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
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()


# Function to create known latent states z
# this puts in NAs everywhere not caught and at first capture
#  we condition on first cpature, and it puts a one in there in the model code, so don't put it in here.
# This is giving some of z so we don't have to estimate it. This z matrix is used in model code, and only the NA's are estimated. Will speed up code to not have to estimate things we know. (It would come up with the right answer becuase of constraints, but would take longer.)
known.state.ms <- function(ms, notseen){
  # notseen: label for 'not seen'
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){ #individual
    m <- min(which(!is.na(state[i,]))) 
    state[i,m] <- NA
  }
  return(state)
}

# Function to create initial values for unknown z
# randomly samples from alive states for every time an animal wasn't caught
ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE) #in current example, 3 states, 3rd is dead
  known.states <- 1:(states-1) # remove unknown state
  v <- which(ch==states) #which ones are not caught?
  ch[-v] <- NA        # all the known states make NA
  ch[v] <- sample(known.states, length(v), replace = TRUE) # unknown sample from known
  return(ch)
}




# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

date()
# Call JAGS from R (BRT 8 min)
ms <- jags(jags.data, inits, parameters, "ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(ms, digits = 3)
date()


par(mfrow = c(3, 2), las = 1)
hist(ms$BUGSoutput$sims.list$mean.phi[,1], col = "gray", main = "", xlab = expression(phi[A]), ylim=c(0,1300))
abline(v = phiA, col = "red")	
hist(ms$BUGSoutput$sims.list$mean.phi[,2], col = "gray", main = "", xlab = expression(phi[B]), ylim=c(0,1300), ylab="")
abline(v = phiB, col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,1], col = "gray", main = "", xlab = expression(psi[AB]), ylim=c(0,1300))
abline(v = psiAB, col="red")
hist(ms$BUGSoutput$sims.list$mean.psi[,2], col = "gray", main = "", xlab = expression(psi[BA]), ylab="", ylim=c(0,1300))
abline(v = psiBA, col="red")
hist(ms$BUGSoutput$sims.list$mean.p[,1], col = "gray", main = "", xlab = expression(p[A]), ylim=c(0,1300))
abline(v = pA, col = "red")
hist(ms$BUGSoutput$sims.list$mean.p[,2], col = "gray", main = "", xlab = expression(p[B]), ylab="", ylim=c(0,1300))
abline(v = pB, col = "red")

######################################################
#################################### alternative model without individual covariates (remove the i in phi and psi arrays)

# Specify model in BUGS language
sink("ms.alternative1.jags")
cat("
    model {
    
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
    }
    for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1)      # Priors for mean state-spec. survival
    mean.psi[u] ~ dunif(0, 1)      # Priors for mean transitions
    mean.p[u] ~ dunif(0, 1)        # Priors for mean state-spec. recapture
    }
    
    # Define state-transition and observation matrices
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    ps[1,t,1] <- phiA[t] * (1-psiAB[t])
    ps[1,t,2] <- phiA[t] * psiAB[t]
    ps[1,t,3] <- 1-phiA[t]
    ps[2,t,1] <- phiB[t] * psiBA[t]
    ps[2,t,2] <- phiB[t] * (1-psiBA[t])
    ps[2,t,3] <- 1-phiB[t]
    ps[3,t,1] <- 0
    ps[3,t,2] <- 0
    ps[3,t,3] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,t,1] <- pA[t]
    po[1,t,2] <- 0
    po[1,t,3] <- 1-pA[t]
    po[2,t,1] <- 0
    po[2,t,2] <- pB[t]
    po[2,t,3] <- 1-pB[t]
    po[3,t,1] <- 0
    po[3,t,2] <- 0
    po[3,t,3] <- 1
    } #t
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()


######################################################
############################################# if don't need temporal effects:
# Specify model in BUGS language
sink("ms.alternative2.jags")
cat("
    model {
    
    # Priors and constraints	
    phiA ~ dunif(0, 1)    # Prior for mean survival in A
    phiB ~ dunif(0, 1)    # Prior for mean survival in B
    psiAB ~ dunif(0, 1)   # Prior for mean movement from A to B
    psiBA ~ dunif(0, 1)   # Prior for mean movement from B to A
    pA ~ dunif(0, 1)      # Prior for mean recapture in A
    pB ~ dunif(0, 1)      # Prior for mean recapture in B
    
    # Define state-transition and observation matrices
    # Define probabilities of state S(t+1) given S(t)
    ps[1,1] <- phiA * (1-psiAB)
    ps[1,2] <- phiA * psiAB
    ps[1,3] <- 1-phiA
    ps[2,1] <- phiB * psiBA
    ps[2,2] <- phiB * (1-psiBA)
    ps[2,3] <- 1-phiB
    ps[3,1] <- 0
    ps[3,2] <- 0
    ps[3,3] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,1] <- pA
    po[1,2] <- 0
    po[1,3] <- 1-pA
    po[2,1] <- 0
    po[2,2] <- pB
    po[2,3] <- 1-pB
    po[3,1] <- 0
    po[3,2] <- 0
    po[3,3] <- 1
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1],])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t],])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# if want to do in parallel, use jags.parallel() instead of jags(). This sends chains to different cores.
# need to pass all the info to jags.data, including functions for inits and known states, number chains, etc