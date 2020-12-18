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

ninds <- c(800,400,400)
ntot <- sum(ninds)
n.months <- 50   
start.date <- sw.temp$date2[20] # Aug 1994
end.date <- start.date + months(n.months-1)
n.webs <- length(ninds)
n.sec.occ <- 3


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
beta.0             = 0.001
beta.male          = 0.005
beta.site          = c(0.01,0.02) 
immig              = 0.01

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
for(w in 1:n.webs){
  obs[[w]]<- data.frame(ID=NA, month=NA, Sec=NA, State=NA)
  
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
                                               State = obs.state))  
        
   } #m
 } #i
  obs[[w]] <- obs[[w]][-1,] # remove the first row of NAs
} #w


# make like SW data
# remove all state = 3
# add indiv covariates like sex, site, web, letter_2, date, session, snv_pos
# dates need to be "1994-08-01", "1994-08-02", "1994-08-03"


####-------------------------------------------------------------------------------###
# Now format the data for each site separately


covariate.data_1 <- # need  cleaned.data, CH.secondary (as list), tags, by.sitetag=TRUE,
  # diversity.data for I.dat  (need MNI though it's not what I used? I guess I'll just use zI)







####-------------------------------------------------------------------------------###
# Now combine the data using the combine data functions


