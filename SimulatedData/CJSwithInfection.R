

#############################################################
## Kery & Schaub (2012), Ch 7 mark-recap, CJS
# Adding infection as a 0 or 1 
#############################################################


##### psi is the prob of remaining uninfected, so (1-psi) is the prob of becoming infected (force of infection)
### so we need to organize our data with 1's for uninfected and 0's for infected 
### here we're assuming that our infection data are correct

# but still need to simulate a z.inf like our z (we don't know the true state when non caught-- except if caught after known infected or non caught between two uninfected captures... so need to fill this in with z.inf matrix somehow)


#### when simulate data, make sure that infection process is embedding with simulating capture histories, because don't want them to contradict each other.

##########################simulate data

# Define parameter values
n.occasions=6						#number of capture occasions
marked=rep(50,n.occasions-1)			#annual number of newly marked indiv

phi=rep(0.65,n.occasions-1)
p=rep(0.4, n.occasions-1)
psi=rep(0.1, n.occasions-1)


#define matrices with survival and recap probs
PHI=matrix(phi,ncol=n.occasions-1,nrow=sum(marked))
P=matrix(p,ncol=n.occasions-1,nrow=sum(marked))
PSI=matrix(psi,ncol=n.occasions-1,nrow=sum(marked))


#define function to simulate a catpure history matrix (CH)
simul.cjs=function(PHI,P,marked){
	n.ocassions=dim(PHI)[2]+1
	CH=matrix(0,ncol=n.occasions,nrow=sum(marked))
	
	#define a vector with the occasion of marking
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

#### now have capture histories, need to say whether infected or uninfected based on psi probabilities

simul.infection=function(CH,PSI){
	n.ocassions=dim(PSI)[2]+1
	
	CH.i=replace(CH,which(CH==0),NA)
	### need to replace 1's with 0's or 1's depending on infection. Need to set initial values- whether first cuaght or not based on proportions that should be infected at that force of infection. Once we have inital values then you should be able to simulate using same arguments for survival using psi instead of phi
	# rough estimate of about 20% infected

	# create vector with occasion of marking (first seen)
#	get.first=function(x) min(which(x!=0))
#	f=apply(CH,1,get.first)
	#each of these first captures is set to uninfected(1) or infected(0) with 80/20% prob
#	for(i in 1:length(f)){
#		CH.i[i,f[i]]=rbinom(1,1,0.8)
		
#	}


#I think it will be easier to simulate true data with infection and survival process, then simulate dataset of capture based on prob of detectio after that to make sure we are simulating infection correctly, when animals weren't caught.	
	

}











#specify model in BUGS language
sink("cjs-c-c.bug")
#cat "					######<--------------------- uncomment 
model{
	
###############Priors and constraints (specifiying the model)
for(i in 1:nind){
	for(t in f[i]:(n.occasions-1)){
		phi[i,t]<-mean.phi	#### here is where could make more complicated by putting in glm
		p[i,t]<-mean.p
		psi[i,t]<-mean.psi
		} #t
	} #i

mean.phi~dunif(0,1)		#prior for mean survival
mean.p~dunif(0,1)		#prior for mean recapture
mean.psi~dunif(0,1)		#prior for probablity of remaining uninfected


#############Likelihood 		## shouldn't have to change this to run diff models
for(i in 1:ind){
	#define latent state at first capture
	z[i,f[i]]=1		# z is true (latent) state alive or dead, know alive at first capture
	for(t in  (f[i]+1):n.occasions){
		
		#state process				# alive or dead
		z[i,t]~dbern(mu1[i,t]) 		#mu1 is probability alive
		mu1[i,t]=phi[i,t-1]*z[i,t-1] # prob alive = survival prob * if it was alive last time (if wasn't alive then multiplying by zero, so this assures that animals that are dead stay dead)

		y.inf[i,t]~dbern(mu2[i,t]) 		#mu1 is probability alive
		mu2[i,t]=psi[i,t-1]*y.inf[i,t-1] # prob uninfected given uninfected at previous time step = psi * if it was uninfected last time (if was infected then multiplying by zero, so this assures that animals that are infected stay infected)
		
		
		####need to put in z.inf (to fill in where didn't capture)
						
		#observation process			# caught or not
		y[i,t]~dbern(mu3[i,t]) 		# mu2 is prob of capture
		mu3[i,t]=p[i,t-1]*z[i,t]	# capture prob= p * if it was alive
		} #t
	} #i
}
#",fill=TRUE)  #####<----------------uncomment this
sink()


