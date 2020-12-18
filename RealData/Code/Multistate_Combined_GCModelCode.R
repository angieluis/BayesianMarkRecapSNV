## Multistate Model specification Combined 3 sites using max model for Grand Canyon


## To do: ----------------------------------------------------------##
# 1:    Error when run, see the bottom. I think I need to make 
#       known state z[i,f] = NA. 
# 2:    I'm assuming can't go from infected to uninfected. Check data.
# 3:    Think about how to model beta. Need a better transform
#         - should vary by site or web
#         - need infected immigration
#         - eventually want to vary by diversity, etc
#         - maybe estimate q
## -----------------------------------------------------------------##


 # specify model
   sink("Multistate_Combinedsites_GCModel.bug")
   
   cat("
     
     model {
       
    # -------------------------------------------------
       # States (S):
       # 1 alive as S
       # 2 alive as I
       # 3 dead

       # Observations (O):  
       # 1 seen as S 
       # 2 seen as I
       # 3 not seen
       # -------------------------------------------------



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
      alpha.inf        ~ dnorm(0, 0.4)T(-10, 10)    # prior for infection
      alpha.inf.male   ~ dnorm(0, 0.4)T(-10, 10)    # prior for interaction with sex (male)

      ## webs
      for(w in 1:(max(web) - 1)) {
        alpha.web[w] ~ dnorm(0, 0.4)T(-10, 10)      # prior for coef on web
        }
      alpha.web.use <- c(0, alpha.web)
      

      ## seasonal covariates
      # don't estimate a beta coefficient for winter (assume 0),
      # then estimate a coefficent for all other seasons
      for(m in 1:3) {
        alpha.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on season
        alpha.ndvi.season[m] ~ dnorm(0, 0.4)T(-10, 10)
        alpha.prcp.season[m] ~ dnorm(0, 0.4)T(-10, 10)
        alpha.tmin.season[m] ~ dnorm(0, 0.4)T(-10, 10)
        alpha.tmax.season[m] ~ dnorm(0, 0.4)T(-10, 10)
       }
      alpha.season.use <- c(0, alpha.season) 
      alpha.ndvi.season.use <- c(0, alpha.ndvi.season)
      alpha.prcp.season.use <- c(0, alpha.prcp.season)
      alpha.tmin.season.use <- c(0, alpha.tmin.season)
      alpha.tmax.season.use <- c(0, alpha.tmax.season)
      
       
      ##### PRIORS FOR PSI #####
      for (s in 1:2){     #for each sex
        sex.beta[s] ~ dunif(0, 1)  # uniform - where the highest beta could be probability 
                                   # of becoming infected=1 at highest MNI (scaled, so 1)
                                # female first, then male
       } 
       
       

      ##### PRIORS FOR RECAPTURE #####
      sigma.0     ~ dnorm(0, 0.4)T(-10, 10)  # prior for intercept on p
      sigma.recap ~ dnorm(0, 0.4)T(-10, 10)  # prior for coeff on recap same occ
      sigma.male  ~ dnorm(0, 0.4)T(-10, 10)  # prior for coef on sex=male
      sigma.inf   ~ dnorm(0, 0.4)T(-10, 10)  # prior for infection
      sigma.inf.male ~ dnorm(0, 0.4)T(-10, 10)  #prior for infection on males (interaction)
      
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
            
            ### Phi for uninfected
            logit(phiS[i, m]) <-            
              
              alpha.0 +
              alpha.male * sex[i] +                 # 0 if female, 1 if male
              alpha.season.use[season[i,m]] +       
              alpha.web.use[web[i]] +
              
              # environmental covariates
              alpha.ndvi * ndvi[i, m] +
              alpha.prcp * prcp[i, m] +
              alpha.tmin * tmin[i, m] +
              alpha.tmax * tmax[i, m] +
              alpha.swe  * swe[i, m] +
              alpha.swe.winter * swe.winter[i, m] +
              alpha.ndvi.season.use[season[i,m]] * ndvi[i, m] +
              alpha.prcp.season.use[season[i,m]] * prcp[i, m] +
              alpha.tmin.season.use[season[i,m]] * tmin[i, m] +
              alpha.tmax.season.use[season[i,m]] * tmax[i, m]

            ### Phi for infected
            logit(phiI[i, m]) <-              
              
              alpha.0 +
              alpha.male * sex[i] +                 # 0 if female, 1 if male
              alpha.season.use[season[i,m]] +       
              alpha.web.use[web[i]] +
       
              # environmental covariates
              alpha.ndvi * ndvi[i, m] +
              alpha.prcp * prcp[i, m] +
              alpha.tmin * tmin[i, m] +
              alpha.tmax * tmax[i, m] +
              alpha.swe  * swe[i, m] +
              alpha.swe.winter * swe.winter[i, m] +
              alpha.ndvi.season.use[season[i,m]] * ndvi[i, m] +
              alpha.prcp.season.use[season[i,m]] * prcp[i, m] +
              alpha.tmin.season.use[season[i,m]] * tmin[i, m] +
              alpha.tmax.season.use[season[i,m]] * tmax[i, m] + 

              # infection terms
              alpha.inf +
              alpha.inf.male * sex[i]
       
            } # m for months
          } # i for individual
        
        
      ##### MODEL FOR PSI #####
        ##### 2 dimensions [indiv, month]

       
      for(i in 1:nind) {

          beta[i] <- sex.beta[sex[i] + 1]         
          
          for(m in 1:(n.months[i] - 1)) {           

            psiSI[i,m] <- beta[i] * I.dat[i, m]   # think about indexing here 
          } #m          
      } #i

        # really should calculate I here based on p (num.caught/p). but I think there will be 
            #estimation problems, so for now, I'll go with using MNI (scaled between (0,1))
       
       # I don't want a logit transform on beta/psi because of the S shape
       # so need to make sure beta is bounded so that beta*max(I.dat)<1.
       # Since I.dat will be scaled with max 1, beta needs to be between 0 and 1.
       # This also doens't account for infected immigration. That means if I[t-1]=0,
       # then no chance of getting infected.  
       # Could also change to mean.beta*I[t-1]/N^q and estimate q to account for 
       # frequency dependent transmission or in between. But then need to worry about the bounds 
       

        ##### MODEL FOR P #####
        ##### 3 dimensions [indiv, month, day]
        for(i in 1:nind) {
          for(m in months.trapped.mat[i, 1:length.months.trapped[i]]) {

            # updated to account for differnt secondary occasions
            for(d in 1:n.sec.occ[i, Prim[i,m]]) {               
              logit(pS[i, m, d]) <-
                sigma.0 +                       # intercept
                sigma.recap * p.or.c[i, m, d] + # if caught before in same session
                sigma.male * sex[i] +           # adj. for males (0 if female)
                sigma.season.use[season[i,m]] + # season factor, where winter=0
                sigma.web.use[web[i]]


              logit(pI[i, m, d]) <-
                sigma.0 +                       
                sigma.recap * p.or.c[i, m, d] +  
                sigma.male * sex[i] +            
                sigma.season.use[season[i,m]] +  
                sigma.web.use[web[i]] +
       
                #### infection terms
                sigma.inf +
                sigma.inf.male * sex[i]
       


              } # d for days
            } # m for month
          } # i for individual
        
      
      ################## Define state-transition and observation matrices
      for (i in 1:nind){  
        # Define probabilities of State(m+1) given State(m)
        for (m in 1:(n.months[i] - 1)){ 
        ps[1,i,m,1] <- phiS[i,m] * (1-psiSI[i,m]) # survival of S to S
        ps[1,i,m,2] <- phiS[i,m] * psiSI[i,m]     # survival and transition from S to I
        ps[1,i,m,3] <- 1-phiS[i,m]                # S to dead
        ps[2,i,m,1] <- 0                          # I to S (can't happen, tho check data)
        ps[2,i,m,2] <- phiI[i,m]                  # survival of I to I
        ps[2,i,m,3] <- 1-phiI[i,m]                # I to dead
        ps[3,i,m,1] <- 0                          # dead to S
        ps[3,i,m,2] <- 0                          # dead to I
        ps[3,i,m,3] <- 1                          # dead stay dead
      } #m
       for (m in months.trapped.mat[i, 1:length.months.trapped[i]]){        
        for(d in 1:n.sec.occ[i, Prim[i,m]]) { 
          # Define probabilities of Observation(m) given State(m)
          # first index is state, last index is observation
          # could potentially include observation as wrong state (false neg or pos)

          po[1,i,m,d,1] <- pS[i,m,d]        # in S and observed as S
          po[1,i,m,d,2] <- 0                # in S and observed as I 
          po[1,i,m,d,3] <- 1-pS[i,m,d]      # in S and not observed 
          po[2,i,m,d,1] <- 0                # in I and observed as S
          po[2,i,m,d,2] <- pI[i,m,d]        # in I and observed as I
          po[2,i,m,d,3] <- 1-pI[i,m,d]      # in I and not observed
          po[3,i,m,d,1] <- 0                # dead and observed as S
          po[3,i,m,d,2] <- 0                # dead and observed as I
          po[3,i,m,d,3] <- 1                # dead and not observed
         } #d
        } #m
       } #i
       

        ##### LIKELIHOOD #####
        
        ##### STATE PROCESS
        for(i in 1:nind) {
          
          # define latent state at first capture 
          # dimensions [individual, month] spans study period
          # z is true (latent) state, know state at first capture
          z[i,f.longmonth[i]] <- f.state[i]    
          
          for(m in (f.longmonth[i] + 1):n.months[i]) {  ### can go backwards if only caught last day?
            z[i, m] ~ dcat(ps[z[i, m-1], i, m-1,])
          } # m (total months spanned)
        } #i
        
        
        ##### OBSERVATION PROCESS
        for(obs in 1:n.obs) {
          y[obs] ~ dcat(po[z[id[obs], longmonth.obs[obs]], id[obs], longmonth.obs[obs], sec[obs],])
        } # obs
        
       ### Compilation error.  z[1,20] is a logical node and cannot be observed
       ### This is the first time that individual was observed. Can I only model from f[i]+1?
       ### But this wasn't a problem with the CJS models....
       
        
        
        } # model
    ", fill = TRUE)
    
    sink()
     
  
## -------------------------------------------------------------------------- ##