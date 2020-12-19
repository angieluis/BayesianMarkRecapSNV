#######################################################################
## Simulate Data for Robust design multistate infection models
## to make sure all the code is working
#######################################################################

#### Right now, missing a few of the complications from the data
#      1) months trapped by site doesn't differ
#      2) all months trapped (longmonth=Prim)
#      3) number secondary occasions always 3
#      4) no p.or.c
#      5) time steps line up (same seasons at same time across sites)
#      6) Simulated data to match obs.dat so when formating data,
#         not starting from the 'dirty data' to to Ch.list 
#         and back to long data. 
######################################################################

# Use sites grandcanyon.e, navajo.1, zuni.1 as a template
# so use covariate data for those sites for the first nmonths

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

ninds <- c(400,200,200)
ntot <- sum(ninds)
n.months <- 50   
start.date <- sw.temp$date2[20] # Aug 1994
end.date <- start.date + months(n.months-1)
n.webs <- length(ninds)
n.sec.occ <- 3
prim.dates <- seq(start.date,end.date,by="month")
sessions <- character()
for(i in 1:length(prim.dates)){
  sessions[i] <- ifelse(month(prim.dates[i])<10 ,paste(year(prim.dates[i]),"0",month(prim.dates[i]),sep=""),paste(year(prim.dates[i]),month(prim.dates[i]),sep=""))
}


### set parameter values
alpha.0            = 1 
alpha.season       = c(1,-0.5,0) 
alpha.web          = c(-0.2,-1.4)
alpha.male         = 0.05 
alpha.swe          = 0.5
alpha.ndvi         = 0.2
alpha.prcp         = 0
alpha.tmin         = -1
alpha.tmax         = 0.5
alpha.ndvi.season  = c(0,0,-0.2)
alpha.prcp.season  = c(0,0.2,0.5)
alpha.tmin.season  = c(0.5,1,0.5)
alpha.tmax.season  = c(0,-1,-0.5)
alpha.swe.winter   = 0
alpha.inf          = -1
alpha.inf.male     = -0.5
sigma.0            = 0.1 #-1.5
#sigma.recap        = 2 
sigma.male         = 0 
sigma.season       = c(0.2,0.2,0.1)
sigma.web          = c(0,0,0)
sigma.inf          = 0.1
sigma.inf.male     = 0.1
beta.0             = 0.00000005
beta.male          = 0.00000005
beta.site          = c(0.00000001,0.00000002) 
immig              = -2

alpha.season.use       = c(0,alpha.season)
alpha.ndvi.season.use  = c(0,alpha.ndvi.season)
alpha.prcp.season.use  = c(0,alpha.prcp.season)
alpha.tmin.season.use  = c(0,alpha.tmin.season)
alpha.tmax.season.use  = c(0,alpha.tmax.season)
sigma.season.use       = c(0,sigma.season)
beta.site.use          = c(0,beta.site)
alpha.web.use          = c(0,alpha.web)
sigma.web.use          = c(0, sigma.web)


## Simulate data for each site separately, then paste them together----------------------#
tempcovdata <- list()
indiv.data <- list()
z <- list()

for(w in 1:n.webs){
  web <- webs[w]
  tempcovdata[[w]] <- sw.temp[which(sw.temp$site.web==web),]
  s.ind <- which(tempcovdata[[w]]$date2==start.date) 
  tempcovdata[[w]] <- tempcovdata[[w]][s.ind:(s.ind+(n.months-1)),]

  n <- ninds[w]
  
  # first half female (0), second half male (1)
  indiv.data[[w]] <- data.frame(ID = 1:n, site.web = web, 
                                sex = c(rep(0,n/2), rep(1,n/2)))
  sex <- indiv.data[[w]]$sex
  
  # Model params as functions of covariates ------------------------#

  ############ Survival
  phiS <- matrix(NA,n,n.months-1)
  phiI <- matrix(NA,n,n.months-1)
  for(i in 1:n) {
    for(m in 1:(n.months - 1)) {            
      
      ### Phi for uninfected
      phiS[i, m] <-  rev.logit(          
        
        alpha.0 +
        alpha.male * sex[i] +                 # 0 if female, 1 if male
        alpha.season.use[tempcovdata[[w]]$season[m]] +       
        alpha.web.use[w] +
        
        # environmental covariates
        alpha.ndvi * tempcovdata[[w]]$ndvi[m] +
        alpha.prcp * tempcovdata[[w]]$prcp[m] +
        alpha.tmin * tempcovdata[[w]]$tmin[m] +
        alpha.tmax * tempcovdata[[w]]$tmax[m] +
        alpha.swe  * tempcovdata[[w]]$swe[m] +
        alpha.swe.winter * tempcovdata[[w]]$swewinter[m] +
        alpha.ndvi.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$ndvi[m] +
        alpha.prcp.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$prcp[m] +
        alpha.tmin.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$tmin[m] +
        alpha.tmax.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$tmax[m] )
      
      ### Phi for infected
      phiI[i, m] <-   rev.logit(           
        
        alpha.0 +
          alpha.male * sex[i] +                 # 0 if female, 1 if male
          alpha.season.use[tempcovdata[[w]]$season[m]] +       
          alpha.web.use[w] +
          
          # environmental covariates
          alpha.ndvi * tempcovdata[[w]]$ndvi[m] +
          alpha.prcp * tempcovdata[[w]]$prcp[m] +
          alpha.tmin * tempcovdata[[w]]$tmin[m] +
          alpha.tmax * tempcovdata[[w]]$tmax[m] +
          alpha.swe  * tempcovdata[[w]]$swe[m] +
          alpha.swe.winter * tempcovdata[[w]]$swewinter[m] +
          alpha.ndvi.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$ndvi[m] +
          alpha.prcp.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$prcp[m] +
          alpha.tmin.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$tmin[m] +
          alpha.tmax.season.use[tempcovdata[[w]]$season[m]] * tempcovdata[[w]]$tmax[m] +
      
        # infection terms
        alpha.inf +
        alpha.inf.male * sex[i])
      
    } # m for months
  } # i for individual
  
  ## Need to calculate psi during simulation
   
   
  # Simulate Z -----------------------------------------------------#
  z[[w]] <- matrix(0,n,n.months)
  # simulate when each indiv first enters pop
  f <- sample(1:n.months,size=n,replace=TRUE)
  # simulate state at first capture. 98% uninfected, 2% infected
  f.st <- sample(1:2,size=n,replace=TRUE,prob=c(0.98,0.02)) 
  
  # add those alive the first month
  alive.first.month <- which(f==1)
  z[[w]][alive.first.month,1] <- f.st[alive.first.month]
  
  for(m in 2: n.months){
    I.lasttime <- sum(length(which(z[[w]][,m-1]==2)) )
    for(i in 1:n){  
      if(m==f[i]){ 
        z[[w]][i,m] <- f.st[i]
      } 
      if(m>f[i]){
        last.state <- z[[w]][i,m-1]
      
        # if last.state uninfected
        if(last.state==1){
          # calculate prob of becoming infected
          beta <- beta.0 + beta.male * sex[i] + beta.site.use[w]
          psiSI <- rev.logit(beta * I.lasttime + immig)
          
          # stay uninfected (1) with prob phiS*(1-psiSI)
          # become infected (2) with prob phiS*psiSI
          # die (0) with prob 1-phiS
          z[[w]][i,m] <- sample(c(1,2,0),size=1,prob=c(
            phiS[i,m-1] * (1-psiSI),
            phiS[i,m-1] * psiSI,
            1-phiS[i,m-1])
          )
        }
        
        # if last.state infected
        if(last.state==2){
         # survive (stay 2) with prob phiI
         # die (0) with prob 1-phiI
          z[[w]][i,m] <- rbinom(1,1,phiI[i,m-1])*2
        }
        #if last.state dead
        if(last.state==0){
          z[[w]][i,m] <- 0
        }
      } # if m>f 
    } #i
  } #m

} #webs



# Simulate observations in a longdata format -----------------------------------------#

# each site separately
obs <- list()
longdata <- list()
for(w in 1:n.webs){
  obs[[w]]<- data.frame(ID=NA, month=NA, Sec=NA, State=NA, sex=NA)
  
  pS <- matrix(NA,ninds[w],n.months)
  pI <- matrix(NA,ninds[w],n.months)

  for(i in 1:ninds[w]){
    f <- min(which(z[[w]][i,]>0))
    
    for(m in f:n.months){
        pS[i, m] <- rev.logit(
          sigma.0 +                       # intercept
          #sigma.recap * p.or.c[i, m, d] + # if caught before in same session
          sigma.male * indiv.data[[w]]$sex[i] +           # adj. for males (0 if female)
          sigma.season.use[tempcovdata[[w]]$season[m]] + # season factor, where winter=0
          sigma.web.use[w])
     
     
        pI[i, m] <- rev.logit(
          sigma.0 +                       
          #sigma.recap * p.or.c[i, m, d] +  
          sigma.male * indiv.data[[w]]$sex[i] +            
          sigma.season.use[tempcovdata[[w]]$season[m]] +  
          sigma.web.use[w] +
       
          #### infection terms
          sigma.inf +
          sigma.inf.male * indiv.data[[w]]$sex[i])
     
        
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


#Need to use I's from z to get I.dat instead of MNI
Isum <- list()
Nsum <- list()
for(w in 1:n.webs){
  Isum[[w]] <- apply(z[[w]],2,function(x){length(which(x==2))})
  Nsum[[w]] <- apply(z[[w]],2,function(x){length(which(x>0))})
}

# turn it into long diversity data
diversity.longdata <- data.frame(site="1", web="0",long.month=1:n.months,session=sessions, year=year(prim.dates),month=month(prim.dates),Prim=1:n.months,MNI=Isum[[1]])
for(w in 2:n.webs){
  diversity.longdata <- rbind(diversity.longdata,data.frame(site=as.character(w), web="0",long.month=1:n.months,session=sessions, year=year(prim.dates),month=month(prim.dates),Prim=1:n.months,MNI=Isum[[w]]))
}

sim.temp <- sw.temp
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[1],"1.0")
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[2],"2.0")
sim.temp$site.web <- replace(sim.temp$site.web,sim.temp$site.web==webs[3],"3.0")
sim.temp$site <- as.character(sim.temp$site)
sim.temp$site <- replace(sim.temp$site,sim.temp$site=="grandcanyon","1")
sim.temp$site <- replace(sim.temp$site,sim.temp$site=="navajo","2")
sim.temp$site <- replace(sim.temp$site,sim.temp$site=="zuni","3")
sim.temp$web <-"0"





####-------------------------------------------------------------------------------###
# Now format the data for each site separately

for(w in 1:n.webs){
  dirty.data <- longdata[[w]]
  cleaned.data <- longdata[[w]]
    cleaned.data$date <- cleaned.data$date.ymd
    
  ms.CH.secondary <- multisite.MS.capture.history.fun(
    dirty.data = dirty.data,
    cleaned.data = cleaned.data, 
    site.webs = c("1.0"),
    species="PM"
  )  

  session.list <- multisite.session.list.fun(
    dirty.data = dirty.data,
    site.webs  = c("1.0")
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
    site.webs = c("1.0"),
    # list of temporal covariates and their time lags
    cov.list = list(prcp = 0, 
                  ndvi = 0, 
                  temp = 0,  
                  tmin = 0,  
                  tmax = 0, 
                  swe  = 0,
                  swewinter = 0,
                  MNI = 0)
  ) ## <------- maybe problem merging diversity data and temp data and/or making covariate matrices

                                        
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
obs.dat$State <- replace(obs.dat$State,obs.dat$State==0,3)

covariate.data <- combine.covariates.fun(covariate.data = list(covariate.data_sim1,covariate.data_sim2,covariate.data_sim3))

monthlyCH <- combine.monthlyCH.fun(monthlyCH = list(monthlyCH_sim1,monthlyCH_sim2,monthlyCH_sim3))

n.sec.occ <- combine.n.sec.occ.func(n.sec.occ = list(n.sec.occ_sim1,n.sec.occ_sim2,n.sec.occ_sim3))

months.trapped.mat <- combine.months.trapped.func(months.trapped.mat = list(months.trapped.mat_sim1,months.trapped.mat_sim2, months.trapped.mat_sim3))

length.months.trapped <- c(length.months.trapped_sim1,length.months.trapped_sim2,length.months.trapped_sim3)


save(obs.dat, covariate.data, monthlyCH, p.or.c, n.sec.occ, months.trapped.mat, length.months.trapped,file="CombinedSimMSData.RData")


