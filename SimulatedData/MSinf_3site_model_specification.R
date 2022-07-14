## Multistate Model specification Combined 3 sites using max model for Grand Canyon


## To do: ----------------------------------------------------------##
## -----------------------------------------------------------------##


 # specify model
   sink("MSinf_3site_model_specification.bug")
   
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
      alpha.0          ~ dnorm(0, 0.4)T(-10, 10)    # prior for survival intercept
      alpha.ndvi       ~ dnorm(0, 0.4)T(-10, 10)    # prior for survival coef on ndvi12
      alpha.inf        ~ dnorm(0, 0.4)T(-10, 10)    # prior for infection on survival
      alpha.inf.male   ~ dnorm(0, 0.4)T(-10, 10)    # prior for infection males on survival


      ## seasonal covariates
      # don't estimate a beta coefficient for winter (assume 0),
      # then estimate a coefficent for all other seasons
      for(m in 1:3) {
        alpha.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on season
       }
      alpha.season.use <- c(0, alpha.season) 

       
      ##### PRIORS FOR PSI #####
      
      beta.0    ~ dnorm(0, 0.4)T(-10, 10)
      beta.male ~ dnorm(0, 0.4)T(-10, 10)
      beta.I    ~ dnorm(0, 1)T(-10, 10)

	  
      ##### PRIORS FOR RECAPTURE #####
      sigma.0     ~ dnorm(0, 0.4)T(-10, 10)  # prior for intercept on p
      sigma.inf   ~ dnorm(0, 0.4)T(-10, 10)  # prior for infection on p
 

        
        ##### MODEL FOR PHI #####
        for(i in 1:nind) {
          for(m in 1:(n.months[i] - 1)) {              # n.months now varies by individual (site)
            # phi has  2 dimensions [indiv, and months]
            
            ### Phi for uninfected
            logit(phiS[i, m]) <-            
              
              alpha.0 +
              alpha.season.use[season[i,m]] +       

              # environmental covariates
              alpha.ndvi * ndvi[i, m] 

            ### Phi for infected
            logit(phiI[i, m]) <-              
              
              alpha.0 +
              alpha.season.use[season[i,m]] +       

              # environmental covariates
              alpha.ndvi * ndvi[i, m] +

              # infection terms
              alpha.inf +
              alpha.inf.male * sex[i]
       
            } # m for months
          } # i for individual
        
        
      ##### MODEL FOR PSI #####
        ##### 2 dimensions [indiv, month]

       
      for(i in 1:nind) {
          
          for(m in 1:(n.months[i] - 1)) {           

            logit(psiSI[i, m]) <- beta.0 + beta.male * sex[i] + beta.I * I.dat[i, m]  
            
            # Should I.dat be [m-1]? No. Psi[1] is the probability of transition over interval 1->2. 
            # So probability of becoming infected by month 2 (Psi[1]) is determined by number of infecteds at month 1 (I.dat[1]).
								 
          } #m          
      } #i

       
       

        ##### MODEL FOR P #####
        ##### 3 dimensions [indiv, month, day]
        for(i in 1:nind) {
          for(m in months.trapped.mat[i, 1:length.months.trapped[i]]) {

            # accounts for different secondary occasions
            for(d in 1:n.sec.occ[i, Prim[i,m]]) {               
              logit(pS[i, m, d]) <-
                sigma.0 


              logit(pI[i, m, d]) <-
                sigma.0 +                       
                sigma.inf 
       


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
        
        
        
        
        } # model
    ", fill = TRUE)
    
    sink()
     
  
## -------------------------------------------------------------------------- ##