##### Simulate NDVI and estimate using weekly time step


library(R2jags)
library(tidyr)
library(dplyr)
setwd("~/Documents/JAGS/BayesianMarkRecapSNV/SimulatedData")
#setwd("~/BayesianMarkRecapSNV/SimulatedData")

source("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/RobustCJSfunctions.R")
#source("~/BayesianMarkRecapSNV/RealData/RobustCJSfunctions.R")

load("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/Zunitemporalcovariates.RData")
#load("Zunitemporalcovariates.RData")


##########################simulate data where survival is a (logit) linear function of NDVI

# Define parameter values
n.occasions=120						#number of capture occasions
marked=rep(3,n.occasions-1)			#annual number of newly marked indiv
set.seed(452);NDVI.dat=matrix(rnorm(n.occasions*3,0,1),ncol=n.occasions,nrow=3)		#simulate normalized precipitation data 
alpha0=0.1							#intercept coefficent for how Prcp affects phi
alpha1=1							#slope coefficient for how prcp affects phi

phi1=revlogit(alpha0+alpha1*NDVI.dat[1,])		
phi2=revlogit(alpha0+alpha1*NDVI.dat[2,])		
phi3=revlogit(alpha0+alpha1*NDVI.dat[3,])		
p=rep(0.4, n.occasions)			# p is constant
c=rep(0.45,n.occasions)

#define matrices with survival and recap probs
PHI1=matrix(phi1,ncol=n.occasions,nrow=sum(marked),byrow=TRUE)
PHI2=matrix(phi2,ncol=n.occasions,nrow=sum(marked),byrow=TRUE)
PHI3=matrix(phi3,ncol=n.occasions,nrow=sum(marked),byrow=TRUE)
P=matrix(p,ncol=n.occasions,nrow=sum(marked))
C=matrix(c,ncol=n.occasions,nrow=sum(marked))


#define function to simulate a catpure history matrix (CH)
simul.cjs.rb <- function(PHI, P, C, marked, n.prim.occasions, n.sec.occasions){
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


sim.data1=simul.cjs.rb(PHI1, P, C, marked, n.prim.occasions=n.occasions, n.sec.occasions=rep(3,n.occasions))
sim.data2=simul.cjs.rb(PHI2, P, C, marked, n.prim.occasions=n.occasions, n.sec.occasions=rep(3,n.occasions))
sim.data3=simul.cjs.rb(PHI3, P, C, marked, n.prim.occasions=n.occasions, n.sec.occasions=rep(3,n.occasions))

CH.secondary <- list()
for(i in 1:n.occasions){
  CH.secondary[[i]] <- rbind(sim.data1$observed.month.list[[i]],sim.data2$observed.month.list[[i]],sim.data3$observed.month.list[[i]])
}

##################################################
 
CH.primary <- primary.ch.fun(CH.secondary)
temporal.covariates <- Zuni.weekly.temporal.covariates
temporal.covariates <- temporal.covariates[1:which(temporal.covariates$Prim==n.occasions),]

individual.covariates <- data.frame(ID=1:dim(CH.primary)[1],web=c(rep(1,dim(sim.data1$observed.month.list[[1]])[1]),rep(2,dim(sim.data2$observed.month.list[[1]])[1]),rep(3,dim(sim.data3$observed.month.list[[1]])[1])))

NDVI1 <- numeric()
NDVI2 <- numeric()
NDVI3 <- numeric()
ints <- diff(which(is.finite(temporal.covariates$Prim)))
for (i in 1:length(ints)){
  NDVI1<-c(NDVI1,rep(NDVI.dat[1,i],ints[i]))
  NDVI2<-c(NDVI2,rep(NDVI.dat[2,i],ints[i]))
  NDVI3<-c(NDVI3,rep(NDVI.dat[3,i],ints[i]))
}
NDVI1[length(NDVI1)+1] <- NDVI.dat[1,dim(NDVI.dat)[2]]
NDVI2[length(NDVI2)+1] <- NDVI.dat[2,dim(NDVI.dat)[2]]
NDVI3[length(NDVI3)+1] <- NDVI.dat[3,dim(NDVI.dat)[2]]
temporal.covariates$NDVI_1 <- NDVI1
temporal.covariates$NDVI_2 <- NDVI2
temporal.covariates$NDVI_3 <- NDVI3


################### Do all the data manipulation in R - create vectors
# and line up all the information and pass to BUGS with bugs.data



#  add the primary occassion to the data
for(i in 1:length(CH.secondary)){
  CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
}


#### define first capture (f)
# z is going to be by week - but capture histories are by month. need f to be in weeks
first.caught <- apply(CH.primary,1,function(x){min(which(x>0))})
individual.covariates$f.week <- temporal.covariates$week[match(first.caught,temporal.covariates$Prim)]




obs.dat <- purrr::map_df(
  CH.secondary, 
  ~ tibble::as_tibble(.x) %>%    # 
    dplyr::mutate(
      ID = 1:n()
    ) %>% 
    tidyr::gather(Sec, State, -Prim, -ID) %>%
    dplyr::select(ID, Prim, Sec, State)
)



#### also need a column that indicates whether that individual
# has been caught before in that primary occasion (do we use p or c?)
p.or.c <- numeric()
for(i in 1:dim(obs.dat)[1]){
  # the times that animal was caught that primary session
  dat <- obs.dat[which(obs.dat$Prim==obs.dat$Prim[i] & obs.dat$ID==obs.dat$ID[i] & obs.dat$State==1),]
  
  if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
    p.or.c[i] <- 0
  }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
    firstcap <- min(as.numeric(dat$Sec))
    p.or.c[i] <- ifelse(firstcap<obs.dat$Sec[i], 1 ,0)   #BUGs doesn't like characters so 0 is p, 1 is c
  }
}
p.or.c <- factor(p.or.c)
obs.dat <- dplyr::mutate(obs.dat, p.or.c)


# don't add a bunch of NA's for weeks not caught
obs.dat.full <- inner_join(obs.dat,temporal.covariates[,c("week","Prim")])
obs.dat.full <- arrange(obs.dat.full,week,ID)





#  Subset observation data to observed bits
obs.dat.full <- dplyr::left_join(obs.dat.full, individual.covariates[,c("ID","f.week")]) %>%
  dplyr::filter(week >= f.week)


#### add individual covariates
obs.dat.full <- dplyr::left_join(obs.dat.full,individual.covariates)





NDVI.w <- rbind(temporal.covariates$NDVI_1, temporal.covariates$NDVI_2, temporal.covariates$NDVI_3)
NDVI <- matrix(NA,ncol=dim(temporal.covariates)[1], nrow=dim(individual.covariates)[1])
for(i in 1:dim(individual.covariates)[1]){
  NDVI[i,] <- NDVI.w[individual.covariates$web[i],]
}

# this creates a weekly capture history to pass into the
# initial values and known state functions
weeklyCH <- weekly.primaryCH.fun(CH.primary,temporal.covariates)

##### Bundle data
bugs.data <- list(
  y = obs.dat.full$State,
  prim = obs.dat.full$Prim,
  sec = obs.dat.full$Sec,
  id = obs.dat.full$ID,
  f = individual.covariates$f.week, 
  p.or.c = obs.dat.full$p.or.c,
  web = individual.covariates$web,
  sex = individual.covariates$sex,
  nind = dplyr::n_distinct(obs.dat.full$ID), 
  n.weeks = max(obs.dat.full$week), 
  #time.int = time.int,
  n.obs = nrow(obs.dat.full),
  weeklyCH = weeklyCH,
  z = known.state.cjs(weeklyCH), 
  cjs.init.z=cjs.init.z, 
  CH.primary=CH.primary,
  covariate.month = temporal.covariates$month,
  week = obs.dat.full$week,
  NDVI = NDVI
) 

#initial values
inits=function(){list(z=cjs.init.z(weeklyCH,f),alpha.month=runif(12,0,1),mean.p=runif(1,0,1),mean.c=runif(1,0,1),alpha.0=runif(1,0,1) ,alpha.NDVI=runif(1,0,1))} 

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c","alpha.0","alpha.month","alpha.NDVI")


date()
Simweekly.rcjs.phi.NDVI.p.c.constant=jags.parallel(data=bugs.data,inits,parameters,"robust_CJS_weekly_phi_NDVI_p_dot_c_dot.bug",n.chains=3,n.thin=6,n.iter=10000,n.burnin=5000)
date() # 

save.image("SimweeklyNDVImodel.RData")

#######################################################################
### estimated alpha values are for weekly survival, but simulated above as monthly survival
# so not directly equivalent (need to simulate weekly to be exact)
phi1.est.w <- revlogit(1.95 + 0.666*NDVI.dat[1,])
phi1.est.m <- phi1.est.w^4 # approx every 4 weeks
plot(phi1,phi1.est.m,ylim=c(0,1),xlim=c(0,1))
abline(0,1)

phi2.est.w <- revlogit(1.95 + 0.666*NDVI.dat[2,])
phi2.est.m <- phi2.est.w^4 # approx every 4 weeks
points(phi2,phi2.est.m)
phi3.est.w <- revlogit(1.95 + 0.666*NDVI.dat[3,])
phi3.est.m <- phi3.est.w^4 # approx every 4 weeks
points(phi3,phi3.est.m)
# pretty close
#I'd day it's working