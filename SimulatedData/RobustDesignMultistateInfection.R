## assumes can't go back from I to S and secondary occasions are closed (can't change states from S to I). So need to check data and make sure that's true.

#############################################################
#
# Robust Design Multistate Infection capture-recapture models
#
##############################################################


# how to organize data? Keep a non-robust design version for y and z?
# where collapsed into primary occasions - this is how will estimate S
# and then separate arrays for each primary session that holds the 
# secondary session info?

# y.secondary[i, d, m] # i=individual, d=day, m=month
# if 60 primary occasions (months) each with 3 secondary occasions (days) and 300 individuals, then 60 matrices that are dimentions 300 by 3.  300,3,60. IF secondary occsions not always same length, then this could be a problem?

# y.primary[i, m] # which is basically summed over days (or ifelse(caught at all,1,0)) so maybe don't need an actual separate array



library(R2jags)
setwd("~/Documents/JAGS")



### need to simulate data- will do later


########################################################################## analyze the data

# Compute vector with occasion of first capture

#### this isn't multi-strata
y.primary <- apply(y.secondary,c(1,3),sum)
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


#########################################################################
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
    # pS: capture probability for S
    # pI: capture probability for I
    # cS: recapture probability (within same primary session) for S
    # cI: recapture probability (within same primary session) for I
    # gammaS: probability of entrance into the population for S (recruited in)
    # gammaI: probability of entrance into the population for I 
    # -------------------------------------------------
    # States (S):
    # 1 not yet alive
    # 2 alive as S
    # 3 alive as I
    # 4 dead
    # Observations (O):  
    # 1 seen as S 
    # 2 seen as I
    # 3 not seen
    # -------------------------------------------------
    


    #### Constraints
    for (t in 1:(n.occasions-1)){
    phiS[t] <- mean.phi[1]
    phiI[t] <- mean.phi[2]
    psiSI[t] <- mean.psi
    pS[t] <- mean.p[1]
    pI[t] <- mean.p[2]
    cS[t] <- mean.c[1]
    cI[t] <- mean.c[2]
    gammaS[t] <- mean.gamma[1]
    gammaI[t] <- mean.gamma[2]
    }
    
    ### Priors
    for (u in 1:2){     #for both states
    mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. capture
    mean.c[u] ~ dunif(0, 1)      # Priors for mean recap within primary session
    mean.gamma[u] ~ dunif(0, 1)  # Priors for prob of entrance
    }
    mean.psi ~ dunif(0, 1)    # Prior for mean transition
    omega ~ dunif(0, 1)       # Prior for including in dataset
    

    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of State(t+1) given State(t)
    for (t in f[i]:(n.occasions-1)){
    ### for proper process model, here say psiSI[t]<-beta*I[t] and define beta in priors instead of psiSI (need a logit or mlogit transform probably? Or is it fine given the priors?)

    ps[1,i,t,1] <- 1 - gammaS[t] - gammaI[t]
    ps[1,i,t,2] <- gammaS[t]
    ps[1,i,t,3] <- gammaI[t]
    ps[1,i,t,4] <- 0
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phiS[t] * (1-psiSI[t])
    ps[2,i,t,3] <- phiS[t] * psiSI[t] 
    ps[2,i,t,4] <- 1-phiS[t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0              # because can't go from I to S (tho check the data)
    ps[3,i,t,3] <- phiI[t]        # ditto
    ps[3,i,t,4] <- 1-phiI[t]
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 1
    
    # Define probabilities of Observation(t) given State(t)
    # could potentially include observation as wrong state (false neg or pos)
    po[1,i,t,1] <- 0
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1
    po[2,i,t,1] <- pS[t]
    po[2,i,t,2] <- 0
    po[2,i,t,3] <- 1-pS[t]
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- pI[t]
    po[3,i,t,3] <- 1-pI[t]
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 1
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
### fill in all known but unobserved state (can't go back to S from I)
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
#for those that are unknown after initial capture, make it alive and in same state as last seen
#  want NAs everywhere we know and 1:f. 
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

# run model 
date()
ms <- jags(bugs.data, inits, parameters, "msSI.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
date() # ~12 minutes

# yay working, and estimates are good.

print(ms, digits = 3)






