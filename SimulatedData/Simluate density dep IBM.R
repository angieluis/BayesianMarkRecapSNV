######################################################################
# simulate individual based model of population dynamics
######################################################################

alpha0 <- 1
alpha1 <- 0.5 

b0 <- 0.8
b1 <- 0.5

time <- 150
N0 <- 50
N.mat <- matrix(0,ncol=time,nrow=N0)
N.mat[,1] <- 1

K=rep(50,time)  # I did a constant K of 50, but should vary as a function of some covariate

for(t in 2:time){
  birth.rate <- alpha0 - (alpha0-alpha1)*sum(N.mat[,t-1])/K[t-1]
  survival.rate <- b0 - (b0-b1) * sum(N.mat[,t-1])/K[t-1]

  # this is the same individual based model we did in class
  N.mat[,t] <- N.mat[,t-1]*rbinom(dim(N.mat)[1],1,survival.rate)
  
  # now we're adding births to it (adding rows to the N.mat for births from existing individuals)
  births <- rpois(sum(N.mat[,t-1]),birth.rate) # number of new rows to add
  babies <- matrix(0,nrow=sum(births),ncol=time) # make more rows
  babies[,t] <- 1 # fill them in with 1's so now there is a new animal
  N.mat <- rbind(N.mat,babies) # add these individuals to existing matrix
}

N <- apply(N.mat,2,sum)
plot.ts(N)  

######################################################################
### Simulate capture histories
######################################################################


# N.mat is the z matrix (true latent state)
# now need to simulate observation on top of this
 
p <- 0.6
CH <- matrix(0,nrow=dim(N.mat)[1],ncol=dim(N.mat)[2])
for (r in 1:dim(N.mat)[1]){
  for(c in 1:dim(N.mat)[2]){
    CH[r,c] <- N.mat[r,c]*rbinom(1,1,p)
  }
}
# delete rows that sum to 0 (animal was never caught)
CH <- CH[-which(apply(CH,1,sum)==0),]



