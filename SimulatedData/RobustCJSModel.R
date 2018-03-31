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
# if 60 primary occasions (months) each with 3 secondary occasions (days) and 300 individuals, then 60 matrices that are dimentions 300 by 3.  300,3,60. IF secondary occsions not always same length, then this could be a problem?

# y.primary[i, m] # which is basically summed over days (or ifelse(caught at all,1,0)) so maybe don't need an actual separate array


########################## need to simulate data where survival is a (logit) linear function of precipitation & p & c are a function of sex

# call full CH with secondary occasions as 2nd dimension and primary as 3rd
# CH.secondary

### code to simulate closed data, going to use it to make robust design secondary occasitons
data.fn <- function(N = 200, T = 3, p = 0.3, c = 0.4){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.eff <- array(NA, dim = N)
  
  # First capture occasion
  yfull[,1] <- rbinom(n = N, size = 1, prob = p)
  
  # Later capture occasions
  for (j in 2:T){
    p.eff <- ifelse(sum(yfull[,1:(j-1)])==0,p,c) #if caught any time previously then c
    # will need to modify for robust design so if caught any time previously in this primary session
    yfull[,j] <- rbinom(n = N, size = 1, prob = p.eff)
  }
  
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, c = c, C = C, T = T, yfull = yfull, yobs = yobs))
}


CH.secondary <- array(0,dim=c(40,3,12))
for(i in 1:12){
  
  #months 1:3
  #indiv 1:10
  for(i in 1:3){
    CH.secondary[,,i] <- rbind(data.fn()$yobs[1:10,],matrix(0,nrow=30,ncol=3))
  } 
  
  #months 4:6
  #indiv 11:20
  for(i in 4:6){
    CH.secondary[,,i] <- rbind(matrix(0,nrow=10,ncol=3),data.fn()$yobs[1:10,],matrix(0,nrow=20,ncol=3))
  } 
  
  #months 7:9
  # indiv 21:30
  for(i in 7:9){
    CH.secondary[,,i] <- rbind(matrix(0,nrow=20,ncol=3),data.fn()$yobs[1:10,],matrix(0,nrow=10,ncol=3))
  } 
  
  #months 10:12
  #indiv 31:40
  for(i in 10:12){
  CH.secondary[,,i] <- rbind(matrix(0,nrow=30,ncol=3),data.fn()$yobs[1:10,])
  } 
}

# [indiv, secondary occasions, primary occasions]

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

#phi <- phi[,-(n.primary.occasions)] #remove last phi because need one less

#############Likelihood 		
for(i in 1:nind){
	#define latent state at first capture # this applies to primary only
  # dimensions [individual, primary session (month)]
	z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
	
	
	
	for(m in  (f[i]+1):n.primary.occasions){
		
	  #state process				# alive or dead
	  z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
		mu1[i, m] <- phi[i, m]*z[i, m-1] # this assures that animals stay dead
		#lukacs lab code has phi[i,m] book has phi[i,m-1]. which one is right? prob doesn't matter. changes indexing so just need to keep track of it.
	
	
	  #observation process			# caught or not
	  #first secondary occasion within a primary occastion:
	  y[i, 1, m] ~ dbern(p.eff[i, 1, m])
	  p.eff[i, 1, m] <- z[i, m] * p[i, 1, m]   #### problem here
	  
	  #loop over rest of secondary occasions per primary
	  for(d in 2:n.secondary.occasions){
		  y[i, d, m] ~ dbern(p.eff[i, d, m]) 		# p.eff is prob of capture
		  p.eff[i, d, m] <- z[i, m] * ifelse(sum(y[i, 1:(d-1), m])==0, p[i, d, m], c[i, d, m])	# capture prob= p * if it was alive
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

### we shouldn't give initial values for those elements of z whose value is specified in the data, they get an NA
#function to create matrix of initial values for latent state z
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
# 20 seconds

#sumarize posteriors
print(robust.cjs,digits=3) #does ok 

# estimating p too high. everything else ok.

traceplot(robust.cjs) 




