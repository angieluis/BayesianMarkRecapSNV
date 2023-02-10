#############################################################################
## Simulate Data and explore the model for density-dependent IPM
##  Robust design model, where environmental drivers (NDVI)
## affect K and that drives survival and birth rates. 
## This version is formatted as array for Jolly-Seber models in
## 'DD_JSMS_model_specification.R'
#############################################################################

setwd("~/Documents/JAGS/BayesianMarkRecapSNV")

######################################################
## Set parameter Values for simulation
######################################################

## Mortality parameters
m0 <- -2.944439 #logit(0.05)  #survival is 0.95 at N=0
me <- -0.8472979 #logit(0.3) # survival at equilibrium is 0.7
m.male <- 0.2 # additional mortality for males


# Recruitment parameters
b0 <- 0.6931472 #log(2) # birth rate at N=0 is 2
# be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
# birth rate at equilibrium is 0.3 = death rate

# K parameters
k.0 <- 3.1        # K is exp(k.0) at mean ndvi in winter 
k.season <- c(0.6,-0.4,-0.3) # spring, summer, fall, (winter is intercept)
k.ndvi <- 0.35

# detection parameters
sigma.0            = 0.1 



#############################################################################
############## Explore relationships
#############################################################################

NK<-seq(0,2,length=50) # N/K, so equilibrium is at 1.

rev.logit <- function(x) {
  exp(x) / (1 + exp(x))
}

ddbirthfun <- function(NK, b0, be){
  exp(b0 + (be-b0)*NK)
}

ddmortfun <- function(NK,m0, me){
  rev.logit(m0 + (me-m0)*NK)
}




birth.rate <- exp(b0 + (be-b0)*NK)

mortality.rate.female <- rev.logit(m0 + (me-m0)*NK)
mortality.rate.male <- rev.logit(m0 + (me-m0)*NK + m.male)


plot(NK,birth.rate,type="l",col="slateblue",ylim=c(0,1.5),
     ylab="Rates",xlab="N/K",main="Density-Dependent Rates")
lines(NK,mortality.rate.female,col="tomato")
lines(NK,mortality.rate.male,col="tomato",lty=2)
abline(v=1,lty=3)
text(1,1,"equilibrium")
legend("topright",c("recruitment rate","mortality rate female",
                    "mortality rate male"),lty=c(1,1,2),
       col=c("slateblue","tomato","tomato"),bty="n",cex=0.8)


### Make K time varying as a function of season and ndvi --------------------#

k.season.use <- c(0, k.season)  # winter, spring, summer, fall

K <- exp(k.0  + k.ndvi*seq(-2,2,by=0.1))
plot(seq(-2,2,by=0.1),K,type="l",col="blue",main="carrying capacity by season and ndvi",xlab="ndvi") 
K2 <- exp(k.0  + k.ndvi*seq(-2,2,by=0.1) + k.season.use[2])
lines(seq(-2,2,by=0.1),K2,col="green")
K3 <- exp(k.0  + k.ndvi*seq(-2,2,by=0.1) + k.season.use[3])
lines(seq(-2,2,by=0.1),K3,col="red")
K4 <- exp(k.0  + k.ndvi*seq(-2,2,by=0.1) + k.season.use[4])
lines(seq(-2,2,by=0.1),K4,col="orange")
legend("topleft",c("winter","spring","summer","fall"),lty=1, col=c("blue","green","red","orange"),bty="n")


##########################################################################
## Simulate data like in "MSinf_3site_SimulateData.R"
# Use sites grandcanyon.e, navajo.1, zuni.1 as a template
# so use covariate data for those sites for the first n.months
##########################################################################


source("RealData/Code/01_sw_data_functions_more.R")
sw.temp <- read.csv("RealData/Data/southwest_covariates_norm.csv")
sw.temp$date2 <- lubridate::dmy(paste("1", sw.temp$date))
sw.temp$site.web <- paste(sw.temp$site,sw.temp$web,sep=".")
# add season covariate 0=winter
sw.temp$season <- ifelse(sw.temp$month == 12 | 
                           sw.temp$month == 1 | 
                           sw.temp$month == 2 | 
                           sw.temp$month == 3, 1, 
                         ifelse(sw.temp$month == 4 | 
                                  sw.temp$month == 5 | 
                                  sw.temp$month == 6, 2, 
                                ifelse(sw.temp$month == 7 | 
                                         sw.temp$month == 8 | 
                                         sw.temp$month == 9, 3, 4)))

webs <- c("grandcanyon.e","navajo.1","zuni.1") 
n.start <- c(20,20,20) ## numebr of individuals to start simulation with at first time step
n.months <- 40   # months to simulate
start.date <- sw.temp$date2[20] # start covariate data at Aug 1994
end.date <- start.date + months(n.months-1)
n.webs <- length(webs)
n.sec.occ <- 3 # secondary occasions
prim.dates <- seq(start.date,end.date,by="month")
sessions <- character()
for(i in 1:length(prim.dates)){
  sessions[i] <- ifelse(month(prim.dates[i])<10 ,paste(year(prim.dates[i]),"0",month(prim.dates[i]),sep=""),paste(year(prim.dates[i]),month(prim.dates[i]),sep=""))
}





## Simulate data for each site separately, then paste them together----------------------#
tempcovdata <- list()
indiv.data <- list()
K <- list()
N <- list()
z <- list() 
Phi.f <- list() # keep track of survival for housekeeping

for(w in 1:n.webs){
  web <- webs[w]
  tempcovdata[[w]] <- sw.temp[which(sw.temp$site.web==web),]
  s.ind <- which(tempcovdata[[w]]$date2==start.date) 
  tempcovdata[[w]] <- tempcovdata[[w]][s.ind:(s.ind+(n.months-1)),]
  # first half female (0), second half male (1)
  
  
  # Simulate Z -----------------------------------------------------#
  z[[w]] <- matrix(NA,n.start[w],n.months)
  zsex <- rep(0:1,n.start[w]/2) # start each web with 1/2 males and 1/2 females
  z[[w]][,1] <- 1 # all n.start are alive at time 1
  N[[w]] <- rep(NA,n.months)
  K[[w]] <- rep(NA,n.months)
  Phi.f[[w]] <-rep(NA,n.months) 
  
  for(m in 2: n.months){
    
    # calculate K for that month based on ndvi and season
    K[[w]][m-1] <- exp(k.0  + k.ndvi*tempcovdata[[w]]$ndvi[m-1] + k.season.use[tempcovdata[[w]]$season[m-1]])
    
    # calculate population size 
    N[[w]][m-1] <- sum(length(which(z[[w]][,m-1]>0)) )
    
    
    for(i in 1:dim(z[[w]])[1]) {
      
      # calculate survival based on sex, and N/K (last month)
      mort <- rev.logit(m0 + (me-m0)* N[[w]][m-1]/ K[[w]][m-1] +
         m.male * zsex[i])
      phi <- 1-mort
      
      last.state <- z[[w]][i,m-1]
      
      #if last.state dead
      if(last.state==0){
        z[[w]][i,m] <- 0  #stay dead
      }
      
      # if last.state alive
      if(last.state==1){
        # stay alive (1) with prob phi
        # die (0) with prob 1-phi
        
        z[[w]][i,m] <- rbinom(1, 1, phi)
      }
      
     } #i
    
    
    # keep track of Phis
    Phi.f[[w]][m-1] <- phi
 
    
    ############## Add individuals (rows) to z based on births
  
    be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
    birth.rate <- exp(b0 + (be-b0) * N[[w]][m-1]/ K[[w]][m-1])
    new.indivs <- rpois(1, birth.rate*sum(length(which(z[[w]][,m]>0)))) 
    z.new <- matrix(rep(c(rep(0, m-1), # before now, the animal was not here
                          1, # now at month m, animal is here
                          rep(NA,n.months-m)), # rest of row is NA
                        new.indivs), byrow=TRUE, nrow=new.indivs, ncol=n.months)
    z[[w]] <- rbind(z[[w]], z.new)
    # assume equal chance of males and females
    zsex <- c(zsex, sample(c(0,1),size=new.indivs,replace=TRUE)) 
    
    
    ############## Save individual data
      indiv.data[[w]] <- data.frame(ID = 1:dim(z[[w]])[1], site.web = web, 
                                    sex = zsex)
      
  } #m
  
} #webs





# take a look at simulated numbers
lapply(z,colSums)


# Simulate observations in an array (could be ragged) -------------------------------# 
# p model
nzrows <- unlist(lapply(z,function(x){dim(x)[1]}))
nzcols <- unlist(lapply(z,function(x){dim(x)[2]}))

p <- array(NA,dim=c(max(nzrows),max(nzcols),n.webs))

for(w in 1:n.webs){
  for(i in 1:nzrows[w]){
    for(m in 1:nzcols[w]){
      p[i, m, w] <- rev.logit(sigma.0 )                       # intercept
    }
  }
}
    

CHlist <- list()
indiv.data.obs <- list()
for(w in 1:n.webs){
  CHweb <- array(NA, dim=c(nzrows[w], nzcols[w], n.sec.occ))
  for(i in 1:nzrows[w]){
    for(m in 1:nzcols[w]){
      
      # simulate 3 days of captures 
      z.state <- z[[w]][i,m]
      
      CHweb[i, m, ] <- rbinom(n.sec.occ, 1, p[i,m,w]) * z.state
    }
    
  }
  # remove rows of array that sum to 0 (animal never caught) and individual data associated with them
  CHlist[[w]] <- CHweb[-which(apply(CHweb,1,sum)==0),,]
  indiv.data.obs[[w]] <- indiv.data[[w]][-which(apply(CHweb,1,sum)==0),]
}

# make it an array
CHlistdims <- lapply(CHlist,dim)
CHlistnrows <- unlist(lapply(CHlistdims,function(x){x[1]}))
CHlistncols <- unlist(lapply(CHlistdims,function(x){x[2]}))
CHarray <- array(NA, dim=c(max(CHlistnrows), max(CHlistncols), n.sec.occ, n.webs))
# also make arrays that are augmented (add fake individuals) and have a dummy occasion (add first occastion when no one has entered yet)
n.aug <- 400
n.add <- n.aug + (max(CHlistnrows)-CHlistnrows)
CHarray.aug.du <- array(NA, dim=c(max(CHlistnrows)+400, max(CHlistncols)+1, n.sec.occ, n.webs))
for(w in 1:n.webs){
  CHarray[1:CHlistnrows[w], 1:CHlistncols[w], , w] <- CHlist[[w]]
  CHarray.aug.du[1:CHlistnrows[w], 2:(CHlistncols[w]+1), , w] <- CHlist[[w]]
  CHarray.aug.du[(CHlistnrows[w]+1):dim(CHarray.aug.du)[1], , , w] <- 0  # make the added individuals have 0 capture histories (need to think about this if later add differences in trapping schedules)
  CHarray.aug.du[ , 1, , w] <- 0 # make the dummy occasion 0
}


#indiv covariates like sex are [i,w] instead of [i]
sex <- matrix(NA, nrow=max(CHlistnrows), ncol= n.webs)
for(w in 1:n.webs){
  sex[1:CHlistnrows[w],w] <- indiv.data.obs[[w]]$sex
}
## covariates like ndvi and season are now [w,m] instead of [i,m]
## for calculating K for dummy time step I need covariate data for a month earlier
## (prob need to update simulation with this too)
ndvi <- matrix(NA, nrow=n.webs , ncol=max(nzcols)+1)
season <- matrix(NA, nrow=n.webs , ncol=max(nzcols)+1)
tempcovdata2 <- list()
for(w in 1:n.webs){
  tempcovdata2[[w]] <- sw.temp[which(sw.temp$site.web==web),]
  s.ind <- which(tempcovdata2[[w]]$date2==start.date) - 1
  tempcovdata2[[w]] <- tempcovdata2[[w]][s.ind:(s.ind+(n.months)),]
  ndvi[w,1:(nzcols[w]+1)] <- tempcovdata2[[w]]$ndvi
  season[w,1:(nzcols[w]+1)] <- tempcovdata2[[w]]$season
}

n.inds <- CHlistnrows

save(CHarray, CHarray.aug.du, n.inds, n.add, sex, ndvi, season,file="SimulatedData/SimDensDepArrayData.RData")

