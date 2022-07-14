#############################################################################
## Simulate Data and explore the model for density-dependent Multistate 
## Infection Robust design model, where environmental drivers (NDVI)
## affect K and that drives survival and birth rates. 


    ## Think about infection - right now have an intercept term that 
    ##   represents density-independent transmission, but it's not infected
    ##   immigration. No new recruits are infected. Maybe get rid of that 
    ##   intercept but include some small chance of new recruits being infected
    ##   (to be estimated) - I worry this will mess with estimation.
#############################################################################

######################################################
## Set parameter Values for simulation
######################################################

## Mortality parameters
m0 <- -2.944439 #logit(0.05)  #survival is 0.95 at N=0
me <- -0.8472979 #logit(0.3) # survival at equilibrium is 0.7
m.male <- 0.2 # additional mortality for males
m.infected <- 0.6

# Recruitment parameters
b0 <- 0.6931472 #log(2) # birth rate at N=0 is 2
# be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))


# K parameters
k.0 <- 3.1        # K is exp(k.0) at mean ndvi in winter 
k.season <- c(0.6,-0.4,-0.3) # spring, summer, fall, (winter is intercept)
k.ndvi <- 0.35

# detection parameters
sigma.0            = 0.1 
sigma.inf          = 0.1

# Infection parameters
#beta.0             = -2 # background infected 
beta.male          = -5
beta.I             = 3
rho                = 0.01 # probability a new recruit is infected (set prior from 0 to 1)  


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
mortality.rate.inffemale <- rev.logit(m0 + (me-m0)*NK + m.infected)
mortality.rate.infmale <- rev.logit(m0 + (me-m0)*NK + m.infected + m.male)


plot(NK,birth.rate,type="l",col="slateblue",ylim=c(0,1.5),
     ylab="Rates",xlab="N/K",main="Density-Dependent Rates")
lines(NK,mortality.rate.female,col="tomato")
lines(NK,mortality.rate.male,col="tomato",lty=2)
lines(NK,mortality.rate.inffemale,col="green")
lines(NK,mortality.rate.infmale,col="green",lty=2)
abline(v=1,lty=3)
text(1,1,"DFE")
legend("topright",c("recruitment rate","mortality rate female uninfected",
                    "mortality rate male uninfected", "mortality rate female 
                    infected","mortality rate male infected"),lty=c(1,1,2,1,2),
       col=c("slateblue","tomato","tomato","green","green"),bty="n",cex=0.8)


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
maxI <- 20 # to standardize I.dat so beta not so small
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
I <- list()
z <- list() 
Phi.uf <- list() # keep track of uninfected female survival for housekeeping

for(w in 1:n.webs){
  web <- webs[w]
  tempcovdata[[w]] <- sw.temp[which(sw.temp$site.web==web),]
  s.ind <- which(tempcovdata[[w]]$date2==start.date) 
  tempcovdata[[w]] <- tempcovdata[[w]][s.ind:(s.ind+(n.months-1)),]
  # first half female (0), second half male (1)
  
  
  # Simulate Z -----------------------------------------------------#
  z[[w]] <- matrix(NA,n.start[w],n.months)
  zsex <- rep(0:1,n.start[w]/2) # start each web with 1/2 males and 1/2 females
  z[[w]][,1] <- c(2,2,rep(1,n.start[w]-2))# start with 1 infected male and 1 infected female # 1 is uninfected, 2 is infected
  N[[w]] <- rep(NA,n.months)
  I[[w]] <- rep(NA,n.months)
  K[[w]] <- rep(NA,n.months)
  Phi.uf[[w]] <-rep(NA,n.months) 
  
  for(m in 2: n.months){
    
    # calculate K for that month based on ndvi and season
    K[[w]][m-1] <- exp(k.0  + k.ndvi*tempcovdata[[w]]$ndvi[m-1] + k.season.use[tempcovdata[[w]]$season[m-1]])
    
    # calculate population size 
    N[[w]][m-1] <- sum(length(which(z[[w]][,m-1]>0)) )
    
    # calculate number infected
    I[[w]][m-1] <- sum(length(which(z[[w]][,m-1]==2)) )
    
    
    for(i in 1:dim(z[[w]])[1]) {
      
      # calculate survival based on sex, infection, and N/K (last month)
      mort <- rev.logit(m0 + (me-m0)* N[[w]][m-1]/ K[[w]][m-1] +
         ifelse(z[[w]][i,m-1]==2, m.infected, 0) +
         m.male * zsex[i])
      phi <- 1-mort
      # calculate prob of transition
        last.state <- z[[w]][i,m-1]
        
        # if last.state uninfected
        if(last.state==1){
          # calculate prob of becoming infected
          psiSI <- rev.logit(#beta.0 + 
            beta.male * zsex[i] + beta.I * I[[w]][m-1]/maxI)
          
          # stay uninfected (1) with prob phi*(1-psiSI) 
          # become infected (2) with prob phi*psiSI
          # die (0) with prob 1-phi
          
          z[[w]][i,m] <- sample(c(1,2,0),size=1,prob=c(
            phi * (1-psiSI),
            phi * psiSI,
            1-phi))
          
          
        }
        
        # if last.state infected
        if(last.state==2){
          # survive (stay 2) with prob phi
          # die (0) with prob 1-phi
          z[[w]][i,m] <- rbinom(1,1,phi)*2
        }
        
        #if last.state dead
        if(last.state==0){
          z[[w]][i,m] <- 0  #stay dead
        }
      
     } #i
    
    
    # keep track of Phis
    Phi.uf[[w]][m-1] <- phi
 
    
    ############## Add individuals (rows) to z based on births
  
    be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
    birth.rate <- exp(b0 + (be-b0) * N[[w]][m-1]/ K[[w]][m-1])
    new.indivs <- rpois(1,birth.rate*sum(length(which(z[[w]][,m]>0)))) 
    z.new <- matrix(rep(c(rep(0,m-1), # before now, the animal was not here
                          sample(c(1,2),1,prob=c(1-rho,rho)), # now at month m, animal is here, could be infected (2) with prob rho
                          rep(NA,n.months-m)),new.indivs),byrow=TRUE,nrow=new.indivs,ncol=n.months)
    z[[w]] <- rbind(z[[w]],z.new)
    # assume equal chance of males and females
    zsex <- c(zsex, sample(c(0,1),size=new.indivs,replace=TRUE)) 
    
    
    ############## Save individual data
      indiv.data[[w]] <- data.frame(ID = 1:dim(z[[w]])[1], site.web = web, 
                                    sex = zsex)
      
  } #m
  
} #webs











# Simulate observations in a longdata format -----------------------------------------#

# each site separately
obs <- list()
ninds <- unlist(lapply(z,function(x){dim(x)[1]})) # number of indiv per web
longdata <- list()
for(w in 1:n.webs){
  obs[[w]]<- data.frame(ID=NA, month=NA, Sec=NA, State=NA, sex=NA)
  
  pS <- matrix(NA,ninds[w],n.months)
  pI <- matrix(NA,ninds[w],n.months)
  
  for(i in 1:ninds[w]){
    f <- min(which(z[[w]][i,]>0))
    
    for(m in f:n.months){
      pS[i, m] <- rev.logit(
        sigma.0 )                       # intercept

      
      pI[i, m] <- rev.logit(
        sigma.0 +                       
          sigma.inf )

      
      # simulate 3 days of captures and paste onto obs
      z.state <- z[[w]][i,m]
      
      p <- ifelse(z.state==1, pS[i,m], ifelse(z.state==2, pI[i,m], 0))
      
      obs.state <- rbinom(n.sec.occ, 1, p) * z.state
      obs.state <- replace(obs.state,obs.state==0, 3) # 3 for not caught 
      
      obs[[w]] <- rbind(obs[[w]], data.frame(ID = rep(i, n.sec.occ), 
                                             month = rep(m, n.sec.occ), 
                                             Sec = 1:n.sec.occ,
                                             State = obs.state,
                                             sex = rep(indiv.data[[w]]$sex[i] ,n.sec.occ)))  
      
    } #m
  } #i
  
  obs[[w]] <- obs[[w]][-1,] # remove the first row of NAs
  longdata[[w]] <- obs[[w]][-which(obs[[w]]$State==3),]
  
  longdata[[w]]$letter_2 <- "pm"
  longdata[[w]]$site <- w
  longdata[[w]]$web <- 0
  longdata[[w]]$tag <- longdata[[w]]$ID
  longdata[[w]]$snv_pos <- longdata[[w]]$State-1
  longdata[[w]]$date.ymd <- prim.dates[longdata[[w]]$month]
  longdata[[w]]$date.ymd <- longdata[[w]]$date.ymd + days(longdata[[w]]$Sec-1)
  longdata[[w]]$date <- paste(month(longdata[[w]]$date.ymd),day(longdata[[w]]$date.ymd), year(longdata[[w]]$date.ymd), sep="/")
  longdata[[w]]$session <- sessions[longdata[[w]]$month]
  
  
} #w

# paste data from all webs together
cleaned.data <- longdata[[1]]
for(w in 2:length(webs)){
  cleaned.data <- rbind(cleaned.data,longdata[[w]])
}
dirty.data <- cleaned.data
cleaned.data$date <- cleaned.data$date.ymd
cleaned.data$site_web <- paste(cleaned.data$site,cleaned.data$web,sep="_")
cleaned.data$tag_site <- paste(cleaned.data$tag,cleaned.data$site,sep="_")
cleaned.data$snv_adj <- cleaned.data$snv_pos


###  make diversity data (not using this because use I from z here ----------#


diversity.list <- list()
site.webs <- sort(unique(cleaned.data$site_web))
site.webs <- sub("_",".",site.webs)
for(i in 1:length(site.webs)){
  site <- unlist(strsplit(site.webs[i],"[.]"))[1]
  web <- unlist(strsplit(site.webs[i],"[.]"))[2]
  diversity.list[[i]] <- diversity.df.function(
    data = cleaned.data,   
    site = site,  
    web = web,  
    sessions = multisite.session.list.fun(dirty.data=dirty.data,
                                          site.webs=site.webs[i])[[1]],  
    interpolate=TRUE, 
    scale=FALSE, 
    include.pm = TRUE)
}
names(diversity.list) <- site.webs

#### Make it long data to match sw.temp.data

diversity.longdata <- diversity.list[[1]]
for(i in 2:length(diversity.list)){
  diversity.longdata <- rbind(diversity.longdata,diversity.list[[i]])
}




## I'm using I from z, comparison: --------------------------------------#
Isum <- list()
Nsum <- list()
for(w in 1:n.webs){
  Isum[[w]] <- apply(z[[w]],2,function(x){length(which(x==2))})
  Nsum[[w]] <- apply(z[[w]],2,function(x){length(which(x>0))})
}

# turn it into long diversity data
diversity.longdata2 <- data.frame(site="1", web="0",long.month=1:n.months,session=sessions, year=year(prim.dates),month=month(prim.dates),Prim=1:n.months,MNI=Isum[[1]])
for(w in 2:n.webs){
  diversity.longdata2 <- rbind(diversity.longdata2,data.frame(site=as.character(w), web="0",long.month=1:n.months,session=sessions, year=year(prim.dates),month=month(prim.dates),Prim=1:n.months,MNI=Isum[[w]]))
}

#MNI and I from z are very close....
plot(Isum[[1]],diversity.list[[1]]$MNI)
abline(0,1)
plot(Isum[[2]],diversity.list[[2]]$MNI)
abline(0,1)
plot(Isum[[3]],diversity.list[[3]]$MNI)
abline(0,1)



## -----------------------------------------------------

sim.temp <- sw.temp
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[1],"1.0")
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[2],"2.0")
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[3],"3.0")
sim.temp$site <- unlist(lapply(strsplit(sim.temp$site.web,"[.]"),function(x){x[1]}))
sim.temp$web <-unlist(lapply(strsplit(sim.temp$site.web,"[.]"),function(x){x[2]}))





####-------------------------------------------------------------------------------###
# Now format the data for each site separately

for(w in 1:n.webs){
  
  ms.CH.secondary <- multisite.MS.capture.history.fun(
    dirty.data = dirty.data,
    cleaned.data = cleaned.data, 
    site.webs = paste(w,".0",sep=""),
    species="PM"
  )  
  
  session.list <- multisite.session.list.fun(
    dirty.data = dirty.data,
    site.webs  = paste(w,".0",sep="")
  )
  
  # makes a primary capture history from secondary occasions
  ms.CH.primary <- primary.MSch.fun(ms.CH.secondary)
  
  # rename secondary occasions from 1 to number of days
  # looking at column names and changing them
  for (m in 1:length(ms.CH.secondary)) {
    colnames(ms.CH.secondary[[m]]) <- 1:dim(ms.CH.secondary[[m]])[2]
  } 
  
  
  # number of secondary occasions - maximum of secondary occasions across webs
  n.sec.occ <- unlist(lapply(ms.CH.secondary, function(x) {
    dim(x)[2]
  }))
  
  
  
  # apply covariate function to normalized data
  covariate.data <- multisite.monthly.covariate.fun(
    cleaned.data = cleaned.data, 
    CH.secondary = ms.CH.secondary,         
    tags = rownames(ms.CH.secondary[[1]]),
    by.sitetag = TRUE, 
    sessions = session.list$all.sessions,
    temporal.data = sim.temp, 
    multistate=TRUE,
    diversity.data = diversity.longdata,
    remove.na=TRUE,
    site.webs = paste(w,".0",sep=""),
    # list of temporal covariates and their time lags
    cov.list = list(prcp = 0, 
                    ndvi = 0, 
                    temp = 0,  
                    tmin = 0,  
                    tmax = 0, 
                    swe  = 0,
                    swewinter = 0,
                    MNI = 0)
  ) 
  
  
  ## Convert to long format ----
  
  obs.dat <- monthly.longdata.CH.fun(ms.CH.secondary, 
                                     covariate.data$temporal.covariates, 
                                     covariate.data$individual.covariates)
  
  
  monthlyCH <- monthly.primary.CH.fun(ms.CH.primary,
                                      covariate.data$temporal.covariates)
  
  
  webmonths <- list()
  # loop through sessions   
  for (i in 1:(length(session.list) - 1)) {
    x <- match(session.list[[i + 1]], 
               covariate.data$temporal.covariates$session)
    webmonths[[i]] <- x[which(is.finite(x))]
  }
  names(webmonths) <- names(session.list)[-1]
  
  
  # matrix for months trapped
  months.trapped.mat <- matrix(NA, 
                               nrow = dim(ms.CH.secondary[[1]])[1], 
                               ncol = max(unlist(lapply(webmonths, length))))
  
  
  # create numeric
  length.months.trapped <- numeric()
  # loop through months trapped
  for (i in 1:dim(months.trapped.mat)[1]) {
    # this is a factor currently
    webnam <- covariate.data$individual.covariates$web[i] 
    webi   <- which(names(webmonths) == paste("web", webnam, sep = "."))
    length.months.trapped[i] <- length(webmonths[[webi]])
    months.trapped.mat[i, 1:length.months.trapped[i]] <- webmonths[[webi]]
  }
  
  # take care of problem of differing number of secondary occasions
  n.sec.occ <- matrix(NA,
                      nrow = dim(ms.CH.primary)[1],
                      ncol = dim(ms.CH.primary)[2])
  
  for(t in 1:dim(n.sec.occ)[2]) {
    n.sec.occ[,t] <- apply(ms.CH.secondary[[t]],
                           1,
                           function(x){length(which(is.finite(x)))})
  }
  
  
  
  assign(paste("ms.CH.secondary_sim", w, sep=""),ms.CH.secondary)
  assign(paste("session.list_sim",w,sep=""),session.list)
  assign(paste("ms.CH.primary_sim",w,sep=""),primary.MSch.fun(ms.CH.secondary))
  assign(paste("n.sec.occ_sim",w,sep=""),n.sec.occ)
  assign(paste("covariate.data_sim",w,sep=""),covariate.data)
  assign(paste("obs.dat_sim",w,sep=""),obs.dat)
  assign(paste("monthlyCH_sim",w,sep=""),monthlyCH)
  assign(paste("webmonths_sim",w,sep=""),webmonths)
  assign(paste("months.trapped.mat_sim",w,sep=""),months.trapped.mat)
  assign(paste("length.months.trapped_sim",w,sep=""),length.months.trapped)
  assign(paste("n.sec.occ_sim",w,sep=""),n.sec.occ)
  
}



####-------------------------------------------------------------------------------###
# Now combine the data using the combine data functions

obs.dat <- combine.obsdat.fun(obs.dats = list(obs.dat_sim1,obs.dat_sim2,obs.dat_sim3), individual.covariates = list(covariate.data_sim1$individual.covariates,covariate.data_sim2$individual.covariates,covariate.data_sim3$individual.covariates))

## remove NA rows of obs.dat (when not trapped that day)
obs.dat <- obs.dat %>%
  dplyr::filter(!is.na(State))

## replace State 0 (not observed) to 3 for multistate infection model
obs.dat$State <- replace(obs.dat$State, obs.dat$State == 0, 3)

covariate.data <- combine.covariates.fun(covariate.data = list(covariate.data_sim1,covariate.data_sim2,covariate.data_sim3))

monthlyCH <- combine.monthlyCH.fun(monthlyCH = list(monthlyCH_sim1,monthlyCH_sim2,monthlyCH_sim3))

n.sec.occ <- combine.n.sec.occ.func(n.sec.occ = list(n.sec.occ_sim1,n.sec.occ_sim2,n.sec.occ_sim3))

months.trapped.mat <- combine.months.trapped.func(months.trapped.mat = list(months.trapped.mat_sim1,months.trapped.mat_sim2, months.trapped.mat_sim3))

length.months.trapped <- c(length.months.trapped_sim1,length.months.trapped_sim2,length.months.trapped_sim3)

# modeled so no c
p.or.c <- array(0, dim=c(dim(months.trapped.mat), 3))

save(obs.dat, covariate.data, monthlyCH, n.sec.occ, p.or.c, months.trapped.mat, length.months.trapped,file="SimulatedData/SimDensDepMSinfData.RData")

