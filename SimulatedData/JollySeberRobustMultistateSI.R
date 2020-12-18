#########################################################################################
#
# Estimation of survival, recruitment and population size using the Jolly-Seber model (multi-state formulation)
## Robust design with data as ragged array
# 
##########################################################################################
library(R2jags)
library(dplyr)

logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}



##################################
#  simulate data

######### Still working on this - need to calculate I for psi=beta*I
##################################

# Define parameter values
n.prim.occasions <- 20						            # number of primary capture occasions
n.sec.occasions <- sample(3:5,n.prim.occasions,replace=TRUE)  # number of secondary occasions, can vary between 3 and 5
N <- 400                                      # Superpopulation size
phiS <- rep(0.8, n.prim.occasions)           #survivial of Susceptible
phiI <- rep(0.7, n.prim.occasions)          #survival of Infected
beta <- rep(0.08, n.prim.occasions)         # transmission rate
#psiSI <- rep(0.3, n.prim.occasions)         # prob of trasition from S to I
gammaS <- rep(0.11, n.prim.occasions)         # Entry probability of S
gammaI <- rep(0.02, n.prim.occasions)   # Entry probability of I
pS <- rep(0.5, n.prim.occasions)               # Capture probabilities
pI <- rep(0.6, n.prim.occasions)               # Capture probabilities



n.states <- 4           # states: not yet entered, S, I, & dead
n.obs <- 3              # observed as S, I or not observed

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state at time t
# Dimension 2: state at time t+1
# Dimension 3: individual
# Dimension 4: time

I.start <- 2 #assuming starting with 2 infected # should probably assign these to individuals rather than just say how many.


PSI.STATE <- array(NA, dim=c(n.states, n.states, N, n.prim.occasions)) # a transition matrix for each individual (i) over time (t)
psiSI <- rep(NA, n.prim.occasions) #
I <- numeric()

PSI.OBS <- array(NA, dim=c(n.states, n.obs, N, n.prim.occasions))

z <- array(0, dim = c(N, n.prim.occasions+1)) # z is actual state (add 1 that will remove later)
y <- list() # y is Ch observed so list of months with secondary occasions (add a dummy occasion of first month with no captures)

for (t in 1:(n.prim.occasions)){
  for (i in 1:N){
    
    
    # 1. State process matrix
    psiSI[t] <- ifelse(t==1, beta[t]*I.start, beta[t]*I[t-1]) 
    PSI.STATE[,,i,t] <- matrix(c(
      1-gammaS[t]-gammaI[t], gammaS[t],            gammaI[t],         0, 
      0,                     phiS[t]*(1-psiSI[t]), phiS[t]*psiSI[t],  1-phiS[t],
      0,                     0,                    phiI[t],           1-phiI[t],
      0,                     0,                    0,                 1       ), 
      nrow = n.states, byrow = TRUE)

    # 2.Observation process matrix
    PSI.OBS[,,i,t] <- matrix(c(
      #seen in S, I, not seen
      0,     0,     1,                  # not entered yet
      pS[t], 0,     1-pS[t],               # actually in S
      0,     pI[t], 1-pI[t],               # actually in I
      0,     0,     1                   # actually dead
    ), nrow = n.states, byrow = TRUE)
  
    
    # STATE PROCESS
    # draw S(t) given S(t-1)
    z[i,t] <- ifelse(t==1, 1, sum(rmultinom(1,1,PSI.STATE[z[i,t-1], ,i, t-1]) *1:4))
    # setting latent state to 1 (not yet entered at 1st occasion, and otherwise go through and figure out the latent state based on state transition matrix)
    
    # OBSERVATION PROCESS 
    # draw O(t) given S(t)
    # z is one longer than y, need to get rid of dummy occasion
    ########## -------- check these subscripts. I'm not sure.
    yt <- matrix(0, nrow=N, ncol=n.sec.occasions[t])
    
    yt[i,] <- colSums(rmultinom(n.sec.occasions[t],1, PSI.OBS[z[i,t], ,i,t-1]) * 1:3) # can't refer to z[t+1] here, because haven't calculated it yet. but can't refer to t-1 either becuase wont work for first time .... really need to get rid of first time step so fill in t+1 instead? 
    y[[t]] <- yt
    
    } #i
  I[t]<- ifelse(t==1,I.start, length(which(z[,t]==3)))  
} #t




  # Remove individuals never captured
  
  
  y2 <- lapply(y,function(x){replace(x,x==3,0)})
  z2 <- replace(z,z==1|z==4,0)
  z2 <- replace(z2,z2==2|z2==3,1)
   
  cap.sum <-  rowSums(matrix(unlist(lapply(y2,rowSums)),nrow=N,ncol=n.prim.occasions,byrow=FALSE))
  never <- which(cap.sum == 0)
  CH.sec <- lapply(y,function(x){x[-never,]})
  Nt <- colSums(z2)    # Actual population size
  
sim <- list(true.state=z[,-1], true.N=Nt[-1], observed.month.list=CH.sec)
  
 

#############################################################
### Analyze data
#############################################################

# secondary capture history as a list of month matrices
CH.secondary <- sim$observed.month.list

# create a primary CH from the secondary capture history:
x <- lapply(CH.secondary,function(x){apply(x,1,min)})
v1 <- unlist(x)
CH.primary <- matrix(v1, nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 

#  Analysis of the JS model as a multistate robust design model
# Add dummy occasion
CH.primary.du <- cbind(rep(3, dim(CH.primary)[1]), CH.primary)
CH.secondary.du <- c(list(matrix(3,nrow=dim(CH.secondary[[1]])[1],ncol=dim(CH.secondary[[1]])[2])), CH.secondary)

# Augment data
nz <- 500
CH.primary.ms <- rbind(CH.primary.du, matrix(3, ncol = dim(CH.primary.du)[2], nrow = nz))
CH.secondary.ms <- lapply(CH.secondary.du,function(x){rbind(x,matrix(3, ncol = dim(x)[2], nrow = nz))})



#####################################################


# The JS model as a multistate Robust Design model
# Specify model in BUGS language
sink("js-rd-ms-inf.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phiS: survival probability of susceptible (S)
    # phiI: survival probability of infected (I)
    # psiSI:  probability of becoming infected (from S to I)
    #     now a function of beta*I (proper process model) 
    # (can't go from I to S)
    # pS: recapture probability as S
    # pI: recapture probability as I      # eventually put in c for S and I
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
    for (t in 1:(n.primary.occasions-1)){
      logit(phiS[t]) <- mean.phi[1]
      logit(phiI[t]) <- mean.phi[2]
      logit(psiSI[t]) <- mean.beta*I[t-1] # proper process model for infection
        # think about how the logit transform makes things an S shape. Is that 
        # what I want?
      logit(pS[t]) <- mean.p[1]
      logit(pI[t]) <- mean.p[2]
      gammaS[t] ~ dunif(0, 1) 
      gammaI[t] ~ dunif(0, 1) 
    }
    
    for (u in 1:2){     #for both states
      mean.phi[u] ~ dnorm(0, 0.4)T(-10,10)    # Priors for mean state-spec. survival
      mean.p[u] ~ dnorm(0, 0.4)T(-10,10)      # Priors for mean state-spec. recapture
    }
    mean.beta ~ dnorm(0, 0.4)T(-10,10)  # Prior # change so can be any non-negative value? This should work tho because beta should be small
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of State(t+1) given State(t)
      for (t in 1:(n.primary.occasions-1)){
        ### for proper process model, here say psiSI[t]<-beta*I[t] and define beta in priors instead of psiSI. Need to calculate I with robust design.
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
        # need to add p or c 
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
    
      y[obs] ~ dcat(po[z[id[obs], prim[obs]], id[obs], prim[obs]-1,]) 
    
    } #obs


    # Calculate derived population parameters
    for (t in 1:(n.primary.occasions-1)){
      qgamma[t] <- 1-gammaS[t]-gammaI[t]  # prob of not entering
    }
    cprob[1] <- gammaS[1]+gammaI[1]   # cummulative prob of entering (either as S or I)
    for (t in 2:(n.primary.occasions-1)){
      cprob[t] <- (gammaS[t]+gammaI[t]) * prod(qgamma[1:(t-1)])
    } #t
    psi <- sum(cprob[])            # Inclusion probability
    for (t in 1:(n.primary.occasions-1)){
      b[t] <- (cprob[t] + 0.001) / (psi + 0.001)       # Entry probability
    } #t
    
    for (i in 1:nind){
      for (t in 2:n.primary.occasions){
        Ss[i,t-1] <- equals(z[i,t], 2) #
      } #t
      for (t in 2:n.primary.occasions){
        Is[i,t-1] <- equals(z[i,t], 3)
      } #t
      for (t in 1:(n.primary.occasions-1)){
        dS[i,t] <- equals(z[i,t]-Ss[i,t],0) 
        dI[i,t] <- equals(z[i,t]-Is[i,t],0)
      } #t   
      aliveS[i] <- sum(Ss[i,])
      aliveI[i] <- sum(Is[i,])
    } #i
    
    for (t in 1:(n.primary.occasions-1)){
      S[t]  <- sum(Ss[,t])       # Actual population size of S
      I[t]  <- sum(Is[,t])       # Actual pop size of I
      N[t]  <- S[t] + I[t]       # Actual total pop size
      BS[t] <- sum(dS[,t])       # Number of S entries
      BI[t] <- sum(dI[,t])       # Number of I entries
      B[t]  <- BS[t] + BI[t]     # total number of entries
      fS[t] <- (BS[t] + 0.001)/(N[t] + 0.001)# per capita recruitment rate of S
      fI[t] <- (BI[t] + 0.001)/(N[t] + 0.001)        # per capita recruitment rate of I
    f[t]  <- (B[t] + 0.001)/(N[t] + 0.001)         # total per capita recruitment rate
    } #t
    for (i in 1:nind){
      w[i] <- 1-equals(aliveS[i],0)-equals(aliveI[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size

  }
    ",fill = TRUE)
sink()


######################

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
known.state.SImsJS <- function(ms=CH.primary.ms, notseen=3){ # ms is multistate capture history
  # notseen: label for 'not seen' #here is 3 
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
  dat <- y[which(y$Prim==y$Prim[i] & y$ID==y$ID[i] & (y$Observation==1|y$Observation==2)),]
  
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
  y =  y$Observation, # Observation of animal (1 seen as S, 2 seen as I, 3 not seen)
  prim = y$Prim, # primary occasion
  # sec = y$Sec, # secondary occasion
  id = y$ID, # individual (1:nind)
  # p.or.c = y$p.or.c, # 0 if not caught before in that session, otherwise 1
  n.obs = nrow(y),
  z = known.state.SImsJS(ms=CH.primary.ms)
)

# Specify initial values
jsmsinf.init <- function(ch=CH.primary.ms, nz){ 
  # ch is primary capture histories after augmentation
  # nz is number of rows added for augmentation

  kn.state <- known.state.SImsJS(ms=ch)
  state <- matrix(2,nrow=dim(ch)[1],ncol=dim(ch)[2]) # default is S (2)
  state <- replace(state,!is.na(kn.state),NA)

  for(i in 1:(dim(state)[1]-nz)){
    f <- min(which(is.na(state[i,])))       # before ever caught
    if(f>1){state[i,1:(f-1)] <- 2}              # tried both 1 and 2 here, still get errors
    if(length(which(kn.state[i,]==3))>0){
      maxI <- max(which(kn.state[i,]==3))
      if(maxI<dim(state)[2] ){
        state[i,(maxI+1):dim(state)[2]] <- 3 # all after caught as I are I (3)
      }
    }
  }
  state[(dim(state)[1]-nz+1):dim(state)[1],] <- 1
  state[,1] <- NA #this is specified in likelihood
  return(state)
}  
  
  

inits <- function(){list(mean.phi = runif(2, 0, 1), 
                         mean.p = runif(2, 0, 1), 
                         mean.psi = runif(1, 0, 1),
                         # mean.c = runif(2, 0, 1), 
                         # put in small values for gammas (because breaks o/w)
                         gammaS = rep(0.01, length(unique(y$Prim))),
                         gammaI = rep(0.01, length(unique(y$Prim))),
                         z = jsmsinf.init(CH.primary.ms, nz))}    


# Parameters monitored
parameters <- c(
  "mean.p", 
  # "mean.c",
  "mean.phi",
  "mean.psi",
  "b", "Nsuper", "N", "B", "f")

# MCMC settings
ni <- 1000
nt <- 3
nb <- 400
nc <- 3


date()
js.rd.ms.inf <- jags(jags.data, inits, parameters, "js-rd-ms-inf.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
date() 
 
#Error in jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
#Error in node z[2,3] 
#Node inconsistent with parents

# or z[373,3]

print(js.rd.ms.inf, digits = 3)

###
