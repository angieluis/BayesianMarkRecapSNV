library(R2jags)
setwd("~/Documents/JAGS")


logit=function(x){
	log(x/(1-x))}
revlogit=function(x){
	exp(x)/(1+exp(x))}




#################################################################
## Kery & Schaub (2012), Ch 7 mark-recap, CJS with time covariate
#################################################################


##########################simulate data where survival is a (logit) linear function of precipitation & p is a function of sex

# Define parameter values
n.occasions=20						#number of capture occasions
marked=rep(20,n.occasions-1)			#annual number of newly marked indiv

Prcp=rnorm(n.occasions-1,0,1)		#simulate normalized precipitation data 
alpha0=0.8							#intercept coefficent for how Prcp affects phi
alpha1=1.2							#slope coefficient for how prcp affects phi
phi=revlogit(alpha0+alpha1*Prcp)		#<------ phi is a function of Prcp

sex=sample(c(rep(1,sum(marked)/2),rep(0,sum(marked)/2))) #randomly assign sex to each indiv with 50/50 prob , 1=male, 0=female
p.male=rep(0.6,n.occasions-1)			# prob of capture for males (not time varying)
p.female=rep(0.4,n.occasions-1)			# prob of capture for females

#gamma0=-0.63
#gamma1=1.02
#p=revlogit(gamma0+gamma1*sex) # then can change the specification below - get rid of loop

#define matrices with survival and recap probs
PHI=matrix(phi,ncol=n.occasions-1,nrow=sum(marked),byrow=TRUE)

P=matrix(NA,ncol=n.occasions-1,nrow=sum(marked))
for(i in 1:sum(marked)){
	if(sex[i]==1){
		P[i,]=p.male}
	if(sex[i]==0){
		P[i,]=p.female
	}
}

#define function to simulate a catpure history matrix (CH)
simul.cjs=function(PHI,P,marked){
	n.ocassions=dim(PHI)[2]+1
	CH=matrix(0,ncol=n.occasions,nrow=sum(marked))
	
	#define a vactor with the occasion of marking
	mark.occ=rep(1:length(marked),marked[1:length(marked)])
	
	#fill the CH Matrix
	for(i in 1:sum(marked)){
		CH[i,mark.occ[i]]=1		#put a 1 at the release occasion
		if(mark.occ[i]==n.occasions) next	#starts next iter of loop if only caught once
		for(t in (mark.occ[i]+1):n.occasions){
			#Bernouli trial: does indiv survive?
			sur=rbinom(1,1,PHI[i,t-1])
			if(sur==0) break				# if dead, move to next indiv
			#Bernouli trial: is indiv recaptured?
			rp=rbinom(1,1,P[i,t-1])
			if(rp==1) CH[i,t]=1
			} #t
		} #i
	return(CH)	
}

CH=simul.cjs(PHI,P,marked)


########################## code up model

# create vector with occasion of marking (first seen)
get.first=function(x) min(which(x!=0))
f=apply(CH,1,get.first)

#specify model in BUGS language
sink("cjs-Prcp-sex.bug")
cat("					######<--------------------- uncomment 
model{
	
###############Priors and constraints (specifiying the model)
for(i in 1:nind){
	for(t in f[i]:(n.occasions-1)){
		logit(phi[i,t])<-alpha0 + alpha1 * Prcp[t]   ### do i want this to be t-1?
		logit(p[i,t])<-gamma0 + gamma1 * sex[i]  ### need to set up indicator sex?
		} #t
	} #i

alpha0~dunif(-5,5)	#prior for alpha0 (intercept coeff for effect of Prcp on phi)
alpha1~dunif(-5,5)	#prior for alpha1 (slope)
gamma0~dunif(-5,5)	# prior for gamma0 (intercept for p) - revlogit(p male?)
gamma1~dunif(-5,5)	# prior for gamma1 (slope for p)

#############Likelihood 		## shouldn't have to change this to run diff models
for(i in 1:nind){
	#define latent state at first capture
	z[i,f[i]]=1		# z is true (latent) state alive or dead, know alive at first capture
	for(t in  (f[i]+1):n.occasions){
		
		#state process				# alive or dead
		z[i,t]~dbern(mu1[i,t]) 		#mu1 is probability alive
		mu1[i,t]=phi[i,t-1]*z[i,t-1] # prob alive = survival prob * if it was alive last time (if wasn't alive then multiplying by zero, so this assures that animals that are dead stay dead)
		
		#observation process			# caught or not
		y[i,t]~dbern(mu2[i,t]) 		# mu2 is prob of capture
		mu2[i,t]=p[i,t-1]*z[i,t]	# capture prob= p * if it was alive
		} #t
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

##### Bundle data
bugs.data=list(y=CH,f=f,nind=dim(CH)[1],n.occasions=dim(CH)[2],z=known.state.cjs(CH),Prcp=Prcp,sex=sex)

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
inits=function(){list(z=cjs.init.z(CH,f),alpha0=runif(1,-5,5),alpha1=runif(1,-5,5),gamma0=runif(1,-5,5),gamma1=runif(1,-5,5))}

#parameters monitored
parameters=c("alpha0","alpha1","gamma0","gamma1")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3



date()
## Call JAGS from R
cjs.Prcp.sex=jags(bugs.data,inits,parameters,"cjs-Prcp-sex.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() #tell how long it ran


#sumarize posteriors
print(cjs.Prcp.sex,digits=3) #does ok 


traceplot(cjs.Prcp.sex) 


#test test test

