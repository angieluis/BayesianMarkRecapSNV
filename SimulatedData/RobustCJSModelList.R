
#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## This does not assume the # of secondary occasions was constant so
# capture histories are in list
#################################################################

library(R2jags)
setwd("~/Documents/JAGS")


logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}


######################################################################################
##########  simulate data 
# simulates a capture history that is a list of matrices each element of the list being a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
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
  return(list(true.state=z,observed=y))	
}

sim.data=simul.cjs.rb(PHI, P, C, marked, n.sec.occasions)

CH.secondary <- sim.data$observed

################
# Now create a primary CH from the secondary capture history:
x <- lapply(CH.secondary,rowSums)
v1 <- unlist(x)
CH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 

n.secondary.occasions <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))

# create a vector of first marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH.primary, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
#rCH.primary <- CH.primary          # Recoded CH
#rCH.primary[rCH.primary==0] <- 3
#rCH.secondary <- CH.secondary          # Recoded CH
#rCH.secondary[rCH.secondary==0] <- 3


### This isn't working. Probably the list notation. How do you specify lists in bugs language?

#specify model in BUGS language
sink("robust_cjs_list.bug")
cat("					######<--------------------- uncomment 
    model{
    
    ###############Priors and constraints
    mean.phi ~ dunif(0, 1) # prior for phi
    mean.p ~ dunif(0, 1)   # prior for p
    mean.c ~ dunif(0, 1)   # prior for c
    
    for(i in 1:nind){
      for(m in f[i]:n.primary.occasions){  ### for p need every time for phi need -1
        # phi has only 2 dimensions [indiv, and primary occasions]
        phi[i,m] <- mean.phi   # could specify covariates here
    
        # p and c have 3 dimensions [[primary]][indiv, secondary]		
        for(d in 1:n.secondary.occasions){
          p[[m]][i, d] <- mean.p  # could specify covariates here
          c[[m]][i, d] <- mean.c  
          } #d for days
        } #m for months
    } #i for individual
    
    #phi <- phi[,-(n.primary.occasions)] #remove last phi because need one less. this didn't work. Only works when phi for every primary occasion. I think this is because not estimating p for first month caught when should. see below. So may need to adjust this when want to estimate N (S and I) from closed secondary sessions.
    
    #############Likelihood 		
    for(i in 1:nind){
      # define latent state at first capture 
      # dimensions [individual, primary session (month)]
      z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
    
      for(m in  (f[i]+1):n.primary.occasions){  # because don't start until time 2, don't estimate p and c for first primary session. This may be a problem when want to estimate N later.
    
        #state process				# alive or dead
        z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
        mu1[i, m] <- phi[i, m]*z[i, m-1] # this assures that animals stay dead
        # Lukacs lab code has phi[i,m]. Book has phi[i,m-1]. Which one is right? 
        # Prob doesn't matter. Changes indexing so just need to keep track of it.
    
    
        #observation process			# caught or not
        #first secondary occasion within a primary occastion:
        y[[m]][i, 1] ~ dbern(p.eff[[m]][i, 1])
        p.eff[[m]][i, 1] <- z[i, m] * p[[m]][i, 1]   
    
        #loop over rest of secondaryp occasions per primary
        for(d in 2:n.secondary.occasions){
          y[[m]][i, d] ~ dbern(p.eff[[m]][i, d]) 		# p.eff is prob of capture
          p.eff[[m]][i, d] <- z[i, m] * ifelse(sum(y[[m]][i, 1:(d-1)])==0, p[[m]][i, d], c[[m]][i, d])	
          # capture prob= (p if not caught previously that session or c if was caught that session) 
          # multiply by if it was alive (so can't capture animal that not alive)
          # think about p and phi and indexing. need p for each month and one less phi
    
        } #d
      } #m
    } #i
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

# [indiv, secondary occasions, primary occasions]
##### Bundle data
bugs.data=list(y=CH.secondary,f=f,nind=dim(CH.secondary[[1]])[1],n.secondary.occasions=n.secondary.occasions,n.primary.occasions=length(CH.secondary),z=known.state.cjs(CH.primary))

###### function to create matrix of initial values for latent state z
# we shouldn't give initial values for those elements of z whose value is specified in the data.They get an NA.
cjs.init.z=function(ch,f){
  for(i in 1:dim(ch)[1]){
    if(sum(ch[i,])==1) next
    n2=max(which(ch[i,]==1))
    ch[i,f[i]:n2]=NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]]=NA
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
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs_list.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() #tell how long it ran
# not working because doesn't understand list format in BUGS language

#sumarize posteriors
print(robust.cjs,digits=3)  


traceplot(robust.cjs) 

