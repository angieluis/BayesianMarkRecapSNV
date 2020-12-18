
#### the model for phi
phi[i, m]) <- revlogit(
  
  # non interaction
    alpha.0 +
    alpha.male * sex[i] +         # 0 if female, 1 if male
    alpha.season.use[season[m]] +
    alpha.web.use[web[i]] +
    
    # environmental covariates
    alpha.tmin * tmin[i, m] +
    alpha.tmax * tmax[i, m] +
    alpha.swe  * swe[i, m] +
    alpha.swe.winter * swe.winter[i, m] +
    
)


################# first set values for estimates, use the mean of the values (maybe later can look at spread)
alpha.0 <- mean(.....)
alpha.male <- 
alpha.season.use <- c(0, , , )
alpha.web.use <- c(0, )
  
# coefficients for environmental covariates
alpha.tmin <- 
alpha.tmax <- 
alpha.swe  <-
alpha.swe.winter <- 
  
# seasonal interactions
alpha.tmin.season.use <- c(0, , , )
alpha.tmax.season.use <- c(0, , , )



#################################################
# Calculate Average seasonal dynamics for Web 1
#################################################
ave.tmin.w1 <-   # calculate an average tmin measure for each season, first one is average winter, etc. So this should be a vector of length 4
ave.tmax.w1 <-  # ditto
ave.swe.w1 <-  # ditto
ave.swe.winter.w1 <-  # ditto
  
phi.female.web1.seasonal <- revlogit(alpha.0 + alpha.male*0 #because female
  + alpha.season.use + alpha.web.use[1] + alpha.tmin*ave.tmin.w1 + alpha.tmax*ave.tmax.w1
  + alpha.swe*ave.swe.w1 + alpha.swe*ave.swe.winter.w1 + alpha.tmin.season.use*ave.tmin.w1 +
    alpha.tmax.season.use*ave.tmax.w1)
# should return a vector of length 4, one for each season

phi.male.web1.seasonal <- revlogit(alpha.0 + alpha.male*1 #because male
 + alpha.season.use + alpha.web.use[1] + alpha.tmin*ave.tmin.w1 + alpha.tmax*ave.tmax.w1
 + alpha.swe*ave.swe.w1 + alpha.swe*ave.swe.winter.w1 + alpha.tmin.season.use*ave.tmin.w1 +
   alpha.tmax.season.use*ave.tmax.w1)

plot.ts(phi.female.web1.seasonal,lwd=2,col="magenta")
lines(phi.male.web1.seasonal,lwd=2,col="blue")

### Could do the same for web 2, but probably not necessary


#######################################################
# Explore the effect of changing Tmin for each season
#######################################################

tmins <- seq(-1,1,length=20)   # range of tmins over which to explore
# assume other environmental covariates are at their means as above

# first, how does tmin affect winter survival for females web 1
phi.female.web1.tmin.winter <- revlogit(alpha.0 + alpha.male*0 #because female
     + alpha.season.use[1] # because winter 
     + alpha.web.use[1] # because web 1
     + alpha.tmin*tmins 
     + alpha.tmax*ave.tmax.w1 
     + alpha.swe*ave.swe.w1 + alpha.swe*ave.swe.winter.w1 
     + alpha.tmin.season.use[1]*tmins 
     + alpha.tmax.season.use[1]*ave.tmax.w1)
# should be vector of length 20 (because that's how many tmins there are)

plot(tmins,phi.female.web1.tmin.winter,ylab="phi",xlab="tmin",main="effect of tmin in winter")


### could do this for all environmental covariates

  
 
