#################################################################
## Robust design CJS w/ diff first capture (p) and recapture (c) probs
## This one does not assume the # of secondary occasions was constant
#################################################################

library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS")


logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}


# Data format: list of capture histories (matrices) for each month
#  each element of the list is a primary occasion (month), and the matrix is i by d (rows are individuals and columns are secondary occasions- days)
# each month contains all animals in dataset (not just those caught that month)

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
  
  # Remove individuals never captured
  cap.sum <- rowSums(matrix(unlist(lapply(y,rowSums)),nrow=sum(marked),ncol=n.prim.occasions,byrow=FALSE))
  never <- which(cap.sum == 0)
  CH.sec <- lapply(y,function(x){x[-never,]})
  Nt <- colSums(z)    # Actual population size
  
  

  return(list(true.state=z, captures.true=y, observed.month.list=CH.sec))	
  
  
}

sim.data=simul.cjs.rb(PHI, P, C, marked, n.sec.occasions)


######################################################################################
# end simulating data 
######################################################################################

CH.secondary <- sim.data$observed.month.list

# create a primary CH from the secondary capture history:
x <- lapply(CH.secondary,rowSums)
v1 <- unlist(x)
CH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 

####################################################specify model in BUGS language
sink("robust_cjs_raggedarray.bug")
cat("
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
      mu1[i,f[i]] <- phi[i,f[i]]
    
      for(m in (f[i]+1):n.primary.occasions){  
        z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
        mu1[i, m] <- phi[i, m] * z[i, m-1] # this assures that animals stay dead
        # Lukacs lab code has phi[i,m]. Book has phi[i,m-1]. Which one is right? 
        # Prob doesn't matter. Changes indexing so just need to keep track of it.
      } # m
    } # i
    
    # OBSERVATION PROCESS 
    for(obs in 1:n.obs){   
    
      y[obs] ~ dbern(z[id[obs], prim[obs]] * ifelse(p.or.c[obs]==0, p[id[obs], sec[obs], prim[obs]],       c[id[obs], sec[obs], prim[obs]]) ) 		# 0 represents p, 1 represents c (if caught before that       session)
   
    } #obs
    
  } #model
",fill=TRUE)  #####<----------------uncomment this
sink()
#########################################################################################



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



################### Do all the data manipulation in R - create vectors
# and line up all the information and pass to BUGS with bugs.data

#  A hack using a loop to add the primary occassion to the data
for(i in 1:length(CH.secondary)){
  CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
}


obs_dat <- purrr::map_df(
  CH.secondary, 
  ~ tibble::as_tibble(.x) %>%    #
    dplyr::mutate(
      ID = 1:n()
    ) %>% 
    tidyr::gather(Sec, State, -Prim, -ID) %>%
    dplyr::select(ID, Prim, Sec, State)
)

first_obs <- obs_dat %>% 
  dplyr::group_by(ID) %>% 
  dplyr::summarise(f = min(Prim[State == 1]))

#  Subset observation data to observed bits
y <- dplyr::left_join(obs_dat, first_obs) %>%
  dplyr::filter(Prim >= f)


#### also need a column that indicates whether that individual
# has been caught before in that primary occasion (do we use p or c?)
p.or.c <- numeric()
for(i in 1:dim(y)[1]){
  # the times that animal was caught that primary session
  dat <- y[which(y$Prim==y$Prim[i] & y$ID==y$ID[i] & y$State==1),]
  
  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
    p.or.c[i] <- 0
  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
    firstcap <- min(as.numeric(dat$Sec))
    p.or.c[i] <- ifelse(firstcap<y$Sec[i], 1 ,0)   #BUGs doesn't like characters so 0 is p, 1 is c
  }
}

y <- data.frame(y,p.or.c)

##### Bundle data
bugs.data <- list(
  y = y$State,
  prim = y$Prim,
  sec = y$Sec,
  id = y$ID,
  f = first_obs$f, 
  p.or.c = y$p.or.c,
  nind = dplyr::n_distinct(y$ID), #n.secondary.occasions=n.secondary.occasions, 
  max.secondary.occasions = max(y$Sec), 
  n.primary.occasions = max(y$Prim), 
  n.obs = nrow(y),
  z = known.state.cjs(CH.primary)
) 


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
inits=function(){list(z=cjs.init.z(CH.primary,bugs.data$f),mean.phi=runif(1,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1))}

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3




date()
## Call JAGS from R #12 minutes
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs_raggedarray.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() #tell how long it ran

# estimating both p and c at 0.4 (p should be 0.3)

#sumarize posteriors
print(robust.cjs,digits=3) 


traceplot(robust.cjs) 



