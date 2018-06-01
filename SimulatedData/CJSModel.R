library(R2jags)
setwd("~/Documents/JAGS")

#############################################################
## Kery & Schaub (2012), Ch 7 mark-recap, CJS
#############################################################


##########################simulate data

# Define parameter values
n.occasions=6						#number of capture occasions
marked=rep(50,n.occasions-1)			#annual number of newly marked indiv

phi=rep(0.65,n.occasions-1)
p=rep(0.4, n.occasions-1)

#define matrices with survival and recap probs
PHI=matrix(phi,ncol=n.occasions-1,nrow=sum(marked))
P=matrix(p,ncol=n.occasions-1,nrow=sum(marked))

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
sink("cjs-c-c.bug")
cat("					######<--------------------- uncomment 
model{
	
###############Priors and constraints (specifiying the model)
for(i in 1:nind){
	for(t in f[i]:(n.occasions-1)){
		phi[i,t]<-mean.phi	#### here is where could make more complicated by putting in glm instead of making it constant (intercept)
		p[i,t]<-mean.p
		} #t
	} #i

mean.phi~dunif(0,1)		#prior for mean survival
mean.p~dunif(0,1)		#prior for mean recapture

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
#  I know we condition on first cpature, but why wouldn't you put a 1 in here because we know they're alive and need to multiply by that for the next step. 
# This is giving some of z so we don't have to estimate it. But how does it know it doesn't need to estimate it? In the likelihood above, all i,t for z are included in the code, not just those that are unknown.
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
bugs.data=list(y=CH,f=f,nind=dim(CH)[1],n.occasions=dim(CH)[2],z=known.state.cjs(CH))

#function to create matrix of initial values for latent state z
### we shouldn't give initial values for those elements of z whose value is specified in the data, they get an NA
# but since we filled in those we know were alive in the known.states.cjs shouldn't that information go in here too? Shouldn't we remvoe the ones that we know- don't need to estimate?
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
inits=function(){list(z=cjs.init.z(CH,f),mean.phi=runif(1,0,1),mean.p=runif(1,0,1))}

#parameters monitored
parameters=c("mean.phi","mean.p")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3


## Call WinBUGS from R
#cjs.c.c=bugs(bugs.data,inits,parameters,"cjs-c-c.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=TRUE,bugs.directory=bugs.dir,working.directory=getwd())


## Call JAGS from R
cjs.c.c=jags(bugs.data,inits,parameters,"cjs-c-c.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)



#sumarize posteriors
print(cjs.c.c,digits=3)

traceplot(cjs.c.c)


#I believe this is all three chains after burnin and thinning:
hist(cjs.c.c$BUGSoutput$sims.list$mean.p)
#same as in matrix form with deviance and mean.phi:
hist(cjs.c.c$BUGSoutput$sims.matrix$mean.p)


#test



