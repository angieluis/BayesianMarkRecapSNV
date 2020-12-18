#########################################################################################
#
# Estimation of survival, recruitment and population size using the Jolly-Seber model (multi-state formulation)
## Robust design with data as ragged array
# 
##########################################################################################
library(R2jags)

logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}



##################################
#  simulate data
##################################

# Define parameter values
n.prim.occasions <- 10						            # number of primary capture occasions
n.sec.occasions <- sample(3:5,n.prim.occasions,replace=TRUE)  # number of secondary occasions, can vary between 3 and 5
N <- 200                                      # Superpopulation size
phi <- rep(0.7, n.prim.occasions-1)           # Survival probabilities
b <- c(0.34, rep(0.11, n.prim.occasions-1))   # Entry probabilities 
p <- rep(0.5, n.prim.occasions)               # Capture probabilities
c <- rep(0.6, n.prim.occasions)               # recapture in same primary session

PHI <- matrix(phi, ncol = n.prim.occasions-1, nrow = N, byrow = T)
P <- matrix(p, ncol = n.prim.occasions, nrow = N, byrow = T)
C <- matrix(c, ncol = n.prim.occasions, nrow = N, byrow = T)

# Function to simulate capture-recapture data under the JS model
simul.js.rb <- function(PHI, P, C, b, N, n.sec.occasions){
  B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
  n.prim.occasions <- dim(PHI)[2] + 1
  z <- array(0, dim = c(N, n.prim.occasions)) # z is actual state
  y <- list() # y is Ch observed so list of months with secondary occasions 
  for(m in 1:n.prim.occasions){
    y[[m]] <- matrix(0,nrow=N, ncol=n.sec.occasions[m])
  } 
  # Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (m in 1:n.prim.occasions){
    ent.occ <- c(ent.occ, rep(m, B[m]))
  }
  
  # Simulating survival
  for (i in 1:N){
    z[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
    if (ent.occ[i] == n.prim.occasions) next
    for (m in (ent.occ[i]+1):n.prim.occasions){ 
      # Bernoulli trial: has individual survived occasion?
      sur <- rbinom(1, 1, PHI[i,m-1])
      ifelse (sur==1, z[i,m] <- 1, break)
    } #m
  } #i

  
  # Simulating capture
  for(m in 1:n.prim.occasions){

    # First day (secondary occasion) that primary occasion- must be p
    y[[m]][,1] <- rbinom(n = N, size = 1, prob = P[, m]) * z[ , m]
    
    T <- n.sec.occasions[m]  
  
    # Later secondary capture occasions - choose p or c
    for (d in 2:T){
      for(i in 1:N){
        p.eff <- ifelse(sum(y[[m]][i, 1:(d-1)])==0, P[i , m] * z[i , m], C[i , m] * z[i , m]) #if caught any time previously that month then c
        y[[m]][i,d] <- rbinom(1, size = 1, prob = p.eff)
    
      } #i
    } #d
  } #m

  
  
  # Remove individuals never captured
  cap.sum <- rowSums(matrix(unlist(lapply(y,rowSums)),nrow=N,ncol=n.prim.occasions,byrow=FALSE))
  never <- which(cap.sum == 0)
  CH.sec <- lapply(y,function(x){x[-never,]})
  Nt <- colSums(z)    # Actual population size
  
  
  
  return(list(true.state=z, captures.true=y, observed.month.list=CH.sec))	
  
}  
  
  
  
# Execute simulation function
sim <- simul.js.rb(PHI, P, C, b, N, n.sec.occasions)


#############################################################
### Analyze data
#############################################################

# secondary capture history as a list of month matrices
CH.secondary <- sim$observed.month.list

# create a primary CH from the secondary capture history:
x <- lapply(CH.secondary,rowSums)
v1 <- unlist(x)
CH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 

#  Analysis of the JS model as a multistate robust design model
# Add dummy occasion
CH.primary.du <- cbind(rep(0, dim(CH.primary)[1]), CH.primary)
CH.secondary.du <- c(list(matrix(0,nrow=dim(CH.secondary[[1]])[1],ncol=dim(CH.secondary[[1]])[2])), CH.secondary)

# Augment data
nz <- 500
CH.primary.ms <- rbind(CH.primary.du, matrix(0, ncol = dim(CH.primary.du)[2], nrow = nz))
CH.secondary.ms <- lapply(CH.secondary.du,function(x){rbind(x,matrix(0, ncol = dim(x)[2], nrow = nz))})

#Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.secondary.ms <- lapply(CH.secondary.ms,function(x){replace(x,x==0,2)})            # Not seen = 2, seen = 1


#####################################################


# The JS model as a multistate Robust Design model
# Specify model in BUGS language
sink("js-rd-ms.jags")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # gamma: removal entry probability
    # p: capture probability
    # c: recapture prob (in same primary session)
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
    mean.phi ~ dnorm(0, 0.4)T(-10,10)    # Prior for mean survival
    mean.p ~ dnorm(0, 0.4)T(-10,10)      # Prior for mean capture
    mean.c ~ dnorm(0, 0.4)T(-10,10)      #  prior for mean recapture # not using now
    
    
    for (m in 1:n.primary.occasions){
        
        gamma[m] ~ dunif(0, 1) # Prior for entry probabilities 
        
                
      for(i in 1:nind){
        # phi has 2 dimensions [indiv, and primary occasions]
        logit(phi[i,m]) <- mean.phi
        

        logit(p[i, m]) <- mean.p  # could specify covariates here
        logit(c[i, m]) <- mean.c  
      } #i for individual 
    } #m for months
    
   

    # Define state-transition and observation matrices 	
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
      for (m in 1:(n.primary.occasions-1)){
        ps[1,i,m,1] <- 1-gamma[m]    # not yet alive
        ps[1,i,m,2] <- gamma[m]       # entering (not yet alive to alive)
        ps[1,i,m,3] <- 0              # can't go from not yet alive to dead
        ps[2,i,m,1] <- 0              # can't go from alive to not yet alive
        ps[2,i,m,2] <- phi[i,m]       # prob of survival (alive to alive)
        ps[2,i,m,3] <- 1-phi[i,m]     # prob of dying (alive to dead)
        ps[3,i,m,1] <- 0              # can't go from dead to not yet alive
        ps[3,i,m,2] <- 0              # can't go from dead to alive
        ps[3,i,m,3] <- 1              # dead, stay dead
    
 ######### made 2 separate po matrices, one for p and another for c 
        # Define probabilities of O(t) given S(t)
        po.p[1,i,m,1] <- 0            # prob of capture if not yet alive =0
        po.p[1,i,m,2] <- 1            # prob of not capture if not yet alive =1
        po.p[2,i,m,1] <- p[i,m]         # prob of capture if alive 
        po.p[2,i,m,2] <- 1-p[i,m]       # prob of not capturing
        po.p[3,i,m,1] <- 0            # prob of capture if dead = 0 
        po.p[3,i,m,2] <- 1            # prob of not capturing if dead =1
 
        # Define probabilities of O(t) given S(t)
        po.c[1,i,m,1] <- 0            # prob of capture if not yet alive =0
        po.c[1,i,m,2] <- 1            # prob of not capture if not yet alive =1
        po.c[2,i,m,1] <- c[i,m]         # prob of capture if alive 
        po.c[2,i,m,2] <- 1-c[i,m]       # prob of not capturing
        po.c[3,i,m,1] <- 0            # prob of capture if dead = 0 
        po.c[3,i,m,2] <- 1            # prob of not capturing if dead =1
      } #m
    } #i
    
    ############### Likelihood 
    # STATE PROCESS
    for (i in 1:nind){
      # Define latent state at first occasion
      z[i,1] <- 1   # Make sure that all individuals are in state 1 at t=1
       # No one has entered yet (state 1) at t=1, because when input data above (Ch.primary.du), add a row before the actual data (where no has entered yet)

      for (m in 2:n.primary.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,m] ~ dcat(ps[z[i,m-1], i, m-1,])
        
      } #m
    } #i

    # OBSERVATION PROCESS 
    # draw O(t) given S(t)
    for(obs in 1:n.obs){   
      if(p.or.c[obs] == 0){
        y[obs] ~ dcat(po.p[z[id[obs], prim[obs]], id[obs], prim[obs]-1,]) 
      }else{
        y[obs] ~ dcat(po.c[z[id[obs], prim[obs]], id[obs], prim[obs]-1,])
      }
    } #obs


    # Calculate derived population parameters
    for (t in 1:(n.primary.occasions-1)){
      qgamma[t] <- 1-gamma[t]
    }
    cprob[1] <- gamma[1]
    for (t in 2:(n.primary.occasions-1)){
      cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)]) # prob of anyone entering now, is prob of available (not yet entered) indiv entering now gamma[t] and prob of not entering at any time in the past, which is product of all previous (1-gamma)
    } #t

    psi <- sum(cprob[])            # Inclusion probability # the prob of entering the pop at any time step (even if never caught) - a proportion will never enter so not part of the dataset (one of the non-existent animals in the augmented dataset)

    for (t in 1:(n.primary.occasions-1)){
      b[t] <- cprob[t] / psi      # Entry probability
    } #t
    
    for (i in 1:nind){
      for (t in 2:n.primary.occasions){
        al[i,t-1] <- equals(z[i,t], 2)             # why is this t-1? because starting with t=2 because t=1 was added? But don't do that for d below.
      } #t
      for (t in 1:(n.primary.occasions-1)){
        d[i,t] <- equals(z[i,t]-al[i,t],0)
      } #t   
      alive[i] <- sum(al[i,])
    } #i
    
    for (t in 1:(n.primary.occasions-1)){
      N[t] <- sum(al[,t])        # Actual population size
      B[t] <- sum(d[,t])         # Number of entries
      f[t] <- B[t]/N[t]         # recruitment rate  <----- I added this
    } #t
    for (i in 1:nind){
      w[i] <- 1-equals(alive[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
    }
    ",fill = TRUE)
sink()


######################

########################################### need to write this!!!
known.state.jsrdms <- function(){
  
}




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
p.or.c <- numeric()
for(i in 1:dim(y)[1]){
  # the times that animal was caught that primary session
  dat <- y[which(y$Prim==y$Prim[i] & y$ID==y$ID[i] & y$Observation==1),]
  
  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
    p.or.c[i] <- 0
  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
    firstcap <- min(as.numeric(dat$Sec))
    p.or.c[i] <- ifelse(firstcap<y$Sec[i], 1 ,0)   # 0 is p, 1 is c
  }
}

y <- data.frame(y,p.or.c)


# Bundle data
jags.data <- list(
  n.primary.occasions = max(y$Prim), 
  # max.secondary.occasions = max(y$Sec), 
  nind = dplyr::n_distinct(y$ID), # number of individuals (rows) including real and augmented
  y =  y$Observation, # Observation of animal (1 seen, 2 not seen)
  prim = y$Prim, # primary occasion
  sec = y$Sec, # secondary occasion
  id = y$ID, # individual (1:nind)
  p.or.c = y$p.or.c, # 0 if not caught before that session, otherwise 1
  n.obs = nrow(y)
  #,z = known.state.jsrdms(CH.primary)
)

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

inits <- function(){list(mean.phi = runif(1, 0, 1), 
                         mean.p = runif(1, 0, 1), 
                         mean.c = runif(1, 0, 1), 
                         z = js.multistate.init(CH.primary.du, nz))}    


# Parameters monitored
parameters <- c(
  "mean.p", 
  "mean.c",
  "mean.phi", 
  "b", "Nsuper", "N", "B", "f")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3


date()
js.rd.ms <- jags(jags.data, inits, parameters, "js-rd-ms.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
date() 



# p.or.c not working ... 
#Error parsing model file:
#syntax error on line 93 near "=="

# doesn't recognize "=="? But it has before!


print(js.rd.ms, digits = 3)

### Do I ignore the first time step in the output because it was added? OR is that already removed? (Seems to be already removed from some of the output)
