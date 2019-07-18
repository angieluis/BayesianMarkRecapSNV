# tried to 'fix' the indexing, but not sure if the derived parameters are right.

# The JS model as a multistate Robust Design model
# Specify model in BUGS language
sink("RDMS_phi_sexinf.bug")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phiS: survival probability of susceptible (S)
    # phiI: survival probability of infected (I)
    # psiSI:  probability of becoming infected (from S to I)
    #     now a function of beta*I (proper process model) 
    # pS: recapture probability as S
    # pI: recapture probability as I      # eventually put in c for S and I
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
    
    # Priors and constraints
    for (t in 1:(n.primary.occasions-1)){
      for(i in 1:nind){
        logit(phiS[i,t]) <- 
          alpha.0 +                       # intercept
          alpha.maleS * sex[i]           # sex
         
      }
       for(i in 1:nind){ # same as S except another coefficent
        logit(phiI[i,t]) <- 
          alpha.0 +                       # intercept
          alpha.maleI * sex[i]            # sex (diff from S)
          alpha.inf                       # infection
      }
      psiSI[i,t] <- beta[i]*I.dat[i,t]  # think about indexing here 
        # really should calculate I here based on p (num.caught/p). but I think there will be estimation problems, so for now, I'll go with using MNI (scaled between (0,1))

      # I don't want t a logit transform on beta/psi because of the S shape
      # so need to make sure beta is bounded so that beta*max(I.dat)<1.
      # Since I.dat will be scaled with max 1, beta needs to be between 0 and 1.
      # This also doens't account for infected immigration. That mean if I[t-1]=0 then no chance of getting infected.  
      # Could also change to mean.beta*I[t-1]/N^q and estimate q to account for 
      # frequency dependent transmission or in between
      # beta/psi should also be a function of individual (e.g., males may have higher)
      
    }
    
    for(t in 1:n.primary.occasions){  
      logit(pS[t]) <- mean.p[1] # assuming 1 p for S and 1 for I (not time-varying, sex, etc)
      logit(pI[t]) <- mean.p[2]
    }


    for (u in 1:2){     #for both states
      mean.p[u] ~ dnorm(0, 0.4)T(-10,10)      # Priors for mean state-spec. recapture
    } ## might want to change this to include recap (c), see other code
    
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)
    alpha.inf ~ dnorm(0, 0.4)T(-10,10)
    alpha.maleS ~ dnorm(0, 0.4)T(-10,10)
    alpha.maleI ~ dnorm(0, 0.4)T(-10,10)
    

    for (s in 1:2){     #for each sex
      sex.beta[s] ~ dunif(0, 1)  # uniform - where the highest beta could be probability of becoming infected=1 at highest MNI (scaled, so 1)
        # female first, then male
    } 
    for(i in 1:nind){
      beta[i] <- sex.beta[sex[i]+1]
    }
    
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of State(t+1) given State(t)
      for (t in 1:(n.primary.occasions-1)){
        ps[1,i,t,1] <- phiS[t] * (1-psiSI[t]) # survival of S to S
        ps[1,i,t,2] <- phiS[t] * psiSI[t]     # survival and transition from S to I
        ps[1,i,t,3] <- 1-phiS[t]              # S to dead
        ps[2,i,t,1] <- 0                      # I to S (can't happen, tho check data)
        ps[2,i,t,2] <- phiI[t]                # survival of I to I
        ps[2,i,t,3] <- 1-phiI[t]              # I to dead
        ps[3,i,t,1] <- 0                      # dead to S
        ps[3,i,t,2] <- 0                      # dead to I
        ps[3,i,t,3] <- 1                      # dead stay dead
      }
        for (t in 1:n.primary.occasions){        
        # Define probabilities of Observation(t) given State(t)
        # first index is state, last index is observation
        # could potentially include observation as wrong state (false neg or pos)
        po[1,i,t,1] <- pS[t]    # in S and observed as S
        po[1,i,t,2] <- 0        # in S and observed as I 
        po[1,i,t,3] <- 1-pS[t]  # in S and not observed 
        po[2,i,t,1] <- 0        # in I and observed as S
        po[2,i,t,2] <- pI[t]    # in I and observed as I
        po[2,i,t,3] <- 1-pI[t]  # in I and not observed
        po[3,i,t,1] <- 0        # dead and observed as S
        po[3,i,t,2] <- 0        # dead and observed as I
        po[3,i,t,3] <- 1        # dead and not observed
      } #t
    } #i
    
    ############### Likelihood 
    # STATE PROCESS
    for (i in 1:nind){
      # Define latent state at first occasion
      z[i,1] <- 1   # Make sure that all individuals are in state 1 at t=1
      # No one has entered yet (state 1) at t=1, because when input data above 
      # (Ch.primary.du), add a row before the actual data (where no has entered yet)
    
      for (m in 2:n.primary.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,m] ~ dcat(ps[z[i,m-1], i, m-1,])
    
      } #m
    } #i
    
    # OBSERVATION PROCESS 
    # draw O(t) given S(t)
    for(obs in 1:n.obs){   
    
      y[obs] ~ dcat(po[z[id[obs], prim[obs]], id[obs], prim[obs],]) 
    
    } #obs
    
    

    }
    ",fill = TRUE)
sink()

