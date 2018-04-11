#################################################################
## Rboust design CJS with time and sex covariate
#################################################################

library(R2jags)
setwd("~/Documents/JAGS")


logit=function(x){
	log(x/(1-x))}
revlogit=function(x){
	exp(x)/(1+exp(x))}


# how to organize data? Keep a non-robust design version for y and z?
# where collapsed into primary occasions - this is how will estimate S
# and then separate arrays for each primary session that holds the 
# secondary session info?

# y.secondary[i, d, m] # i=individual, d=day, m=month
# if 60 primary occasions (months) each with 3 secondary occasions (days) and 300 individuals, then 60 matrices that are dimentions 300 by 3.  300,3,60. IF secondary occasions not always same length, then this could be a problem?

# y.primary[i, m] # which is basically summed over days (or ifelse(caught at all,1,0)) so maybe don't need an actual separate array


######################################################################################
##########  simulate data 
### problem. Animals are being detected that are not alive.
######################################################################################

# Define parameter values
n.prim.occasions <- 20						#number of primary capture occasions
n.sec.occasions <- 3            # number of secondary occasions

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
  y <- array(0, dim = c(sum(marked), n.sec.occasions, n.prim.occasions)) # y is Ch observed so includes secondary occasions
  
  
  #define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  
  #fill the CH Matrix
  for(i in 1:sum(marked)){
    z[i, mark.occ[i]] <- 1		#put a 1 at the release occasion
    
    # for first month caught 
    #Bernouli trial: is indiv captured?
    ########  secondary occasions, d for days
    for(d in 1:n.sec.occasions){
    p.eff <- ifelse(sum(y[i, 1:(d-1), mark.occ[i]])==0, P[i, mark.occ[i]], C[i, mark.occ[i]]) #if caught any time previously in this session then use c instead of p
    y[i, d, mark.occ[i]] <- rbinom(1, 1, prob = p.eff)
     } #d
    # if never caught then randomly pick a secondary occasion for capture
    if(sum(y[i, , mark.occ[i]])==0){y[i, sample(1:n.sec.occasions, 1), mark.occ[i]] <- 1}
    
    if(mark.occ[i] == n.prim.occasions) next	#starts next iter of loop if caught only at last occasion
    
    for(m in (mark.occ[i]+1):n.prim.occasions){ # m is primary occasion (month)
      #p.eff <- array(NA, dim = sum(marked))
      mu1 <- PHI[i, m-1] * z[i, m-1] # this assures that animals stay dead
      z[i, m] <-  rbinom(1, 1, mu1) 		#Bernouli trial: does animal survive
      
      if(mu1==0) break				# if dead, move to next indiv
    
      #Bernouli trial: is indiv captured?
      ########  secondary occasions, d for days
      for(d in 1:n.sec.occasions){
        p.eff <- ifelse(sum(y[i, 1:(d-1), m])==0, P[i, m], C[i, m]) * z[i, m] #if caught any time previously in this session (m) then c, and if not alive, can't be caught
        y[i, d, m] <- rbinom(1, 1, prob = p.eff)
      } #d
    } #m
  } #i
  return(list(true.state=z,observed=y))	
}

sim.data=simul.cjs.rb(PHI, P, C, marked, n.sec.occasions)

CH.secondary <- sim.data$observed
######################################################################################
# end simulating data 
######################################################################################




####################################
CH.primary <- apply(CH.secondary,c(1,3),sum)
CH.primary <- replace(CH.primary,CH.primary>1,1)

# create a vector of first marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH.primary, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in S, 2 = seen alive in I, 3 = not seen
#rCH.primary <- CH.primary          # Recoded CH
#rCH.primary[rCH.primary==0] <- 3
#rCH.secondary <- CH.secondary          # Recoded CH
#rCH.secondary[rCH.secondary==0] <- 3



#specify model in BUGS language
sink("robust_cjs.bug")
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

    # p and c have 3 dimensions [indiv, secondary, primary]		
    for(d in 1:n.secondary.occasions){
      p[i, d, m] <- mean.p  # could specify covariates here
      c[i, d, m] <- mean.c  
	   } #d for days
    } #m for months
	} #i for individual

#phi <- phi[,-(n.primary.occasions)] #remove last phi because need one less. this didn't work. Only works when phi for every primary occasion. I think this is because not estimating p for first month caught when should. see below. So may need to adjust this when want to estimate N (S and I) from closed secondary sessions.

#############Likelihood 		
for(i in 1:nind){
	#define latent state at first capture 
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
	  y[i, 1, m] ~ dbern(p.eff[i, 1, m])
	  p.eff[i, 1, m] <- z[i, m] * p[i, 1, m]   
	  
	  #loop over rest of secondaryp occasions per primary
	  for(d in 2:n.secondary.occasions){
		  y[i, d, m] ~ dbern(p.eff[i, d, m]) 		# p.eff is prob of capture
		  p.eff[i, d, m] <- z[i, m] * ifelse(sum(y[i, 1:(d-1), m])==0, p[i, d, m], c[i, d, m])	
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
bugs.data=list(y=CH.secondary,f=f,nind=dim(CH.secondary)[1],n.secondary.occasions=dim(CH.secondary)[2],n.primary.occasions=dim(CH.secondary)[3],z=known.state.cjs(CH.primary))

###### function to create matrix of initial values for latent state z
# we shouldn't give initial values for those elements of z whose value is specified in the data.
# they get an NA
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
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() #tell how long it ran
# 20 min for 20 time steps and 20 marked each time
# 6-8 minutes for  15 time steps and 20 marked each time

#sumarize posteriors
print(robust.cjs,digits=3) #does ok 


traceplot(robust.cjs) 




