## Model specification Combined 3 sites data using max model for Grand Canyon

### Things changed from Code file 08:---------------------------------------- ##
#1: added back in the non-interaction terms for env covariates e.g. alpha.ndvi (incl fixing typo)
#2: n.months, season, and Prim now vary by individual (site)
#3: removed derived params

#4: didn't simulate animals that were only caught on the last day of trapping
#       at that sites, so (f.longmonth[i] + 1):n.months[i]) couldn't go
#       backwards, and removed from obs.dat



     
## Grand Canyon model specification ----------------------------------------- ##
## -------------------------------------------------------------------------- ##


 # specify Grand Canyon model
   sink("Code/Combinedsites_GCModel.bug")
   
   cat("
     
     model {
       
      ##### PRIORS AND CONSTRAINTS #####
      
      
      ##### PRIORS FOR PHI #####
      alpha.0          ~ dnorm(0, 0.4)T(-10, 10)    # prior for intercept
      alpha.male       ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on sex = m 
      alpha.ndvi       ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on ndvi12
      alpha.prcp       ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on prcp6
      alpha.tmin       ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on tmin0
      alpha.tmax       ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on tmax6
      alpha.swe        ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef on swe0
      alpha.swe.winter ~ dnorm(0, 0.4)T(-10, 10)    # prior for coef swe.winter
      

      ## webs
      for(w in 1:(max(web) - 1)) {
        alpha.web[w] ~ dnorm(0, 0.4)T(-10, 10)      # prior for coef on web
        }
      alpha.web.use <- c(0, alpha.web)
      

      ## season
      # don't estimate a beta coefficient for winter (assume 0),
      # then estimate a coefficent for all other seasons
      for(m in 1:3) {
        alpha.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on season
        }
      alpha.season.use <- c(0, alpha.season) 
      

      ## interactions
      # interaction term for season * ndvi
      for(m in 1:3) {
        alpha.ndvi.season[m] ~ dnorm(0, 0.4)T(-10, 10)
      }
      alpha.ndvi.season.use <- c(0, alpha.ndvi.season)
      
      
      # interaction term for season * prcp
      for(m in 1:3) {
        alpha.prcp.season[m] ~ dnorm(0, 0.4)T(-10, 10)
      }
      alpha.prcp.season.use <- c(0, alpha.prcp.season)
      
      
      # interaction term for season * tmin
      for(m in 1:3) {
        alpha.tmin.season[m] ~ dnorm(0, 0.4)T(-10, 10)
      }
      alpha.tmin.season.use <- c(0, alpha.tmin.season)
      
      
      # interaction term for season * tmax
      for(m in 1:3) {
        alpha.tmax.season[m] ~ dnorm(0, 0.4)T(-10, 10)
        }
      alpha.tmax.season.use <- c(0, alpha.tmax.season)
      
      
      ##### PRIORS FOR RECAPTURE #####
      sigma.0     ~ dnorm(0, 0.4)T(-10, 10)  # prior for intercept on p
      sigma.recap ~ dnorm(0, 0.4)T(-10, 10)  # prior for coeff on recap same occ
      sigma.male  ~ dnorm(0, 0.4)T(-10, 10)  # prior for coef on sex=male
      
      
      # season
      for(se in 1:3) {
        sigma.season[se] ~ dnorm(0, 0.4)T(-10, 10)   # prior for season
        }
      sigma.season.use <- c(0, sigma.season)
      
      
      # web  
      for(w in 1:(max(web) - 1)) {
        sigma.web[w] ~ dnorm(0, 0.4)T(-10, 10)       # prior for coef on web
        }
        sigma.web.use <- c(0, sigma.web)
        
        
        ##### MODEL FOR PHI #####
        for(i in 1:nind) {
          for(m in 1:(n.months[i] - 1)) {              # n.months now varies by individual (site)
            # phi has  2 dimensions [indiv, and months]
            
            logit(phi[i, m]) <-
              
              # non interaction
              alpha.0 +
              alpha.male * sex[i] +                 # 0 if female, 1 if male
              alpha.season.use[season[i,m]] +       # season is now a matrix
              alpha.web.use[web[i]] +
              
              # environmental covariates
              alpha.ndvi * ndvi[i, m] +
              alpha.prcp * prcp[i, m] +
              alpha.tmin * tmin[i, m] +
              alpha.tmax * tmax[i, m] +
              alpha.swe  * swe[i, m] +
              alpha.swe.winter * swe.winter[i, m] +

              
              # interactions
              alpha.ndvi.season.use[season[i,m]] * ndvi[i, m] +
              alpha.prcp.season.use[season[i,m]] * prcp[i, m] +
              alpha.tmin.season.use[season[i,m]] * tmin[i, m] +
              alpha.tmax.season.use[season[i,m]] * tmax[i, m]
            
            } # m for months
          } # i for individual
        
        
        ##### MODEL FOR P #####
        ##### 3 dimensions [indiv, month, day]
        for(i in 1:nind) {
          for(m in months.trapped.mat[i, 1:length.months.trapped[i]]) {

            # updated to account for differnt secondary occasions
            for(d in 1:n.sec.occ[i, Prim[i,m]]) {               
              logit(p[i, m, d]) <-
                sigma.0 +             # intercept
                sigma.recap * p.or.c[i, m, d] + 
                # adjustment for if animal was caught previously in this primary 
                # session (0 if not caught before and 1 if so)
                
                sigma.male * sex[i] +         # adj. for males (0 if female)
                sigma.season.use[season[i,m]] + # season factor, where winter=0
                sigma.web.use[web[i]]
              
              } # d for days
            } # m for month
          } # i for individual
        
        
        ##### LIKELIHOOD #####
        
        ##### STATE PROCESS
        for(i in 1:nind) {
          
          # define latent state at first capture 
          # dimensions [individual, month] spans study period
          # z is true (latent) state alive or dead, know alive at first capture
          z[i,f.longmonth[i]] <- 1
          
          for(m in (f.longmonth[i] + 1):n.months[i]) {  ### can go backwards if only caught last day
            mu1[i, m] <- phi[i, m - 1] * z[i, m - 1]
            z[i, m] ~ dbern(mu1[i, m])
            } # m (total months spanned)
          } #i
        
        
        ##### OBSERVATION PROCESS
        for(obs in 1:n.obs) {
          y[obs] ~ dbern(z[id[obs], 
                           longmonth.obs[obs]] * p[id[obs], 
                                                   longmonth.obs[obs], sec[obs]])
          } # obs
        
      
        
        
        } # model
    ", fill = TRUE)
    
    sink()
     
  
## -------------------------------------------------------------------------- ##