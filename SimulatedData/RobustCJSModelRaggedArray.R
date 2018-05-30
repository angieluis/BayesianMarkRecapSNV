#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## This one does not assume the # of secondary occasions was constant
#################################################################

library(R2jags)
setwd("~/Documents/JAGS")


logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}


# Data format: Data frame where
# each row is an observation - each time an individual was caught, like the rodent data
# e.g.  individual, ID1, was caught in the third month on days 2 and 3:
# $individual  $primary   $secondary
#   ID1           3           2
#   ID1           3           3
######################################################################################
##########  simulate data 
######################################################################################

# Define parameter values
n.prim.occasions <- 20						#number of primary capture occasions
n.sec.occasions <- sample(3:5,n.prim.occasions,replace=TRUE)  # number of secondary occasions, can vary between 3 and 5

marked <- rep(20, n.prim.occasions-1)			#annual number of newly marked indiv

phi <- rep(0.65,n.prim.occasions-1)
p <- rep(0.3, n.prim.occasions)
c <- rep(0.4, n.prim.occasions)

#define matrices with survival and recap probs
PHI <- matrix(phi, ncol=n.prim.occasions-1, nrow=sum(marked))
P <- matrix(p, ncol=n.prim.occasions, nrow=sum(marked))
C <- matrix(c, ncol=n.prim.occasions, nrow=sum(marked))

#define function to simulate a catpure history matrix (CH)
simul.cjs.rb <- function(PHI, P, C, marked, n.sec.occasions){
  n.prim.ocassions <- dim(PHI)[2]+1
  z <- array(0, dim = c(sum(marked), n.prim.occasions)) # z is actual state
  y <- list() # y is Ch observed so includes secondary occasions, needs to be a list if number of secondary occasions differs across months
  for(m in 1:n.prim.occasions){
    y[[m]] <- matrix(0,nrow=sum(marked), ncol=n.sec.occasions[m])
  } 
  
  
  #define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  
  #fill the CH Matrix
  for(i in 1:sum(marked)){
    z[i, mark.occ[i]] <- 1		#put a 1 at the release occasion
    
    # for first month caught 
    #Bernouli trial: is indiv captured?
    ########  secondary occasions, d for days
    for(d in 1:n.sec.occasions[mark.occ[i]]){
      p.eff <- ifelse(sum(y[[ mark.occ[i]]][i,1:(d-1)])==0, P[i, mark.occ[i]], C[i, mark.occ[i]]) #if caught any time previously in this session then use c instead of p
      
      y[[mark.occ[i]]][i, d] <- rbinom(1, 1, prob = p.eff)
    } #d
    # if never caught then randomly pick a secondary occasion for capture
    if(sum(y[[mark.occ[i]]][i, ])==0){y[[mark.occ[i]]][i, sample(1:n.sec.occasions[mark.occ[i]], 1)] <- 1}
    
    if(mark.occ[i] == n.prim.occasions) next	#starts next iter of loop if caught only at last occasion
    
    for(m in (mark.occ[i]+1):n.prim.occasions){ # m is primary occasion (month)
      mu1 <- PHI[i, m-1] * z[i, m-1] # this assures that animals stay dead
      z[i, m] <-  rbinom(1, 1, mu1) 		#Bernouli trial: does animal survive
      
      if(mu1==0) break				# if dead, move to next indiv
      
      #Bernouli trial: is indiv captured?
      ########  secondary occasions, d for days
      for(d in 1:n.sec.occasions[m]){
        p.eff <- ifelse(sum(y[[m]][i, 1:(d-1)])==0, P[i, m], C[i, m]) * z[i, m] #if caught any time previously in this session (m) then c, and if not alive, can't be caught
        y[[m]][i, d] <- rbinom(1, 1, prob = p.eff)
      } #d
    } #m
  } #i
  
  
  observations <- data.frame(individual=NA, primary=NA, secondary=NA)
  for(i in 1:sum(marked)){
    second<-numeric()
    prim <- numeric()
    ind <- numeric()
    for(m in 1:length(y)){
      second <- which(y[[m]][i,]==1)
      prim <- rep(m, length(second))
      ind <- rep(i, length(second))
      if(length(second)>0){
        observations <- rbind(observations,data.frame(individual=ind,primary=prim,secondary=second))
      }
    }
  }
  observations <- observations[-1,]
  
  return(list(true.state=z,observed.month.list=y,observed.data=observations))	
  
  
}

sim.data=simul.cjs.rb(PHI, P, C, marked, n.sec.occasions)

observed.data <- sim.data$observed.data

######################################################################################
# end simulating data 
######################################################################################





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
  CH.primary[observed.data[i,1], observed.data[i,2]] <- 1
}

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
#rCH.primary <- CH.primary          # Recoded CH
#rCH.primary[rCH.primary==0] <- 3
#rCH.secondary <- CH.secondary          # Recoded CH
#rCH.secondary[rCH.secondary==0] <- 3



# create a vector of first marking
f <- numeric()
for(i in 1:nind){
  f[i] <- min(observed.data$primary[which(observed.data$individual==i)])
}



#specify model in BUGS language
sink("robust_cjs_raggedarray.bug")
cat("					######<--------------------- uncomment 
model{
  
  ###############Priors and constraints
  mean.phi ~ dnorm(0, 0.4)T(-10,10)     # prior for mean survival
  mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
  mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
  
  for(i in 1:nind){
    for(m in f[i]:n.primary.occasions){  
      
      # phi has only 2 dimensions [indiv, and primary occasions]
      logit(phi[i,m]) <- mean.phi   # could specify covariates here
      
      # p and c have 3 dimensions [indiv, secondary, primary]		
      # will use the largest secondary occasion to determine dimensions and they won't all be used in the data / estimation
      for(d in 1:max.secondary.occasions){
        logit(p[i, d, m]) <- mean.p  # could specify covariates here
        logit(c[i, d, m]) <- mean.c  
      } #d for days
    } #m for months
  } #i for individual
  
  
  #############Likelihood 		
  # STATE PROCESS
  for(i in 1:nind){
    # define latent state at first capture 
    # dimensions [individual, primary session (month)]
    z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
    
    for(m in  (f[i]+1):n.primary.occasions){  
      z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
      mu1[i, m] <- phi[i, m] * z[i, m-1] # this assures that animals stay dead
      # Lukacs lab code has phi[i,m]. Book has phi[i,m-1]. Which one is right? 
      # Prob doesn't matter. Changes indexing so just need to keep track of it.
    } # m
  } # i
  
  # OBSERVATION PROCESS 
  for(obs in 1:n.obs){   
    idays <- y$secondary[which(y$individual==y$individual[obs] & y$primary==y$primary[obs])] # days that individual was caught that month
      
    p.eff <- z[y$individual[obs], y$primary[obs]] * # if it was caught before this day make c, otherwise p
        ifelse(any(idays < y$secondary[obs]), c[y$individual[obs], y$secondary[obs], y$primary[obs]], p[y$individual[obs], y$secondary[obs], y$primary[obs]])	
      
    y[obs,] ~ dbern(p.eff) 		# p.eff is prob of capture
    # think about p and phi and indexing. need p for each month and one less phi

  } #obs
  
}
    ",fill=TRUE)  #####<----------------uncomment this
sink()




#function to create matrix with info about known latent state z
known.state.cjs=function(ch){
  state=ch
  for(i in 1:dim(ch)[1]){
    n1=min(which(ch[i,]==1))
    n2=max(which(ch[i,]==1))
    state[i,n1:n2]=1
    state[i,n1]=NA			#only filling in those that were 0s but we know were alive because caught before and after
  }
  state[state==0]=NA
  return(state)
}


##### Bundle data
bugs.data=list(y=observed.data, f=f, nind=nind, n.secondary.occasions=n.secondary.occasions, max.secondary.occasions=max.secondary.occasions, n.primary.occasions=n.primary.occasions, z=known.state.cjs(CH.primary)) 

###### function to create matrix of initial values for latent state z
# we shouldn't give initial values for those elements of z whose value is specified in the data.
# they get an NA
cjs.init.z=function(ch,f){
  for(i in 1:dim(ch)[1]){
    if(sum(ch[i,])==1) next
    n2=max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)	
}

#initial values
inits=function(){list(z=cjs.init.z(CH.primary,f),mean.phi=runif(1,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1))}

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3



date()
## Call JAGS from R
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs_raggedarray.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() #tell how long it ran


#sumarize posteriors
print(robust.cjs,digits=3) 


traceplot(robust.cjs) 




