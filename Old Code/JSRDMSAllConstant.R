# tried to 'fix' the indexing, but not sure if the derived parameters are right.

# The JS model as a multistate Robust Design model
# Specify model in BUGS language
sink("js-rd-ms-inf.jags")
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
    # gammaS: prob of entry into S class 
    # gammaI: prob of entry into I class (infected immigration)
    # -------------------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive as S
    # 3 alive as I
    # 4 dead
    # Observations (O):  
    # 1 seen as S 
    # 2 seen as I
    # 3 not seen
    # -------------------------------------------------
    
    # Priors and constraints
    for (t in 1:(n.primary.occasions-1)){
      logit(phiS[i,t]) <- mean.phi[1]
      logit(phiI[i,t]) <- mean.phi[2]
      logit(psiSI[i,t]) <- mean.beta*I[t-1] # proper process model for infection
      # think about how the logit transform makes things an S shape. Is that 
      # what I want? PRob not. Need a diff transform. Ask Paul.
      # This also doens't account for infected immigration. Assumes all infection
      # comes locally. 
      # Could also change to mean.beta*I[t-1]/N^q and estimate q to account for 
      # frequency dependent transmission or in between
      # beta should also be a function of individual (e.g., males may have higher)
      
      gammaS[t] ~ dunif(0, 1) # prob of entrance could also vary by i?
      gammaI[t] ~ dunif(0, 1) 
    }
    
    for(t in 1:n.primary.occasions){
      logit(pS[t]) <- mean.p[1]
      logit(pI[t]) <- mean.p[2]
    }

    for (u in 1:2){     #for both states
      mean.phi[u] ~ dnorm(0, 0.4)T(-10,10)    # Priors for mean state-spec. survival
      mean.p[u] ~ dnorm(0, 0.4)T(-10,10)      # Priors for mean state-spec. recapture
    }
    
    mean.beta ~ dnorm(0, 0.4)T(-10,10)  # Prior # change so can be any non-negative value? This should work tho because beta should be small
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of State(t+1) given State(t)
      for (t in 1:(n.primary.occasions-1)){
        ps[1,i,t,1] <- 1-gammaS[t]-gammaI[t]  # still not yet entered
        ps[1,i,t,2] <- gammaS[t]              # just entered as S
        ps[1,i,t,3] <- gammaI[t]              # just entered as I
        ps[1,i,t,4] <- 0                      # not yet entered to dead
        ps[2,i,t,1] <- 0                      # S to not yet entered
        ps[2,i,t,2] <- phiS[t] * (1-psiSI[t]) # survival of S to S
        ps[2,i,t,3] <- phiS[t] * psiSI[t]     # survival and transition from S to I
        ps[2,i,t,4] <- 1-phiS[t]              # S to dead
        ps[3,i,t,1] <- 0                      # I to not yet entered
        ps[3,i,t,2] <- 0                      # I to S (can't happen, tho check the data)
        ps[3,i,t,3] <- phiI[t]                # survival of I to I
        ps[3,i,t,4] <- 1-phiI[t]              # I to dead
        ps[4,i,t,1] <- 0                      # dead to not yet entered
        ps[4,i,t,2] <- 0                      # dead to S
        ps[4,i,t,3] <- 0                      # dead to I
        ps[4,i,t,4] <- 1                      # dead stay dead
      }
        for (t in 1:n.primary.occasions){        
        # Define probabilities of Observation(t) given State(t)
        # first index is state, last index is observation
        # need to add p or c 
        # could potentially include observation as wrong state (false neg or pos)
        po[1,i,t,1] <- 0        # not yet entered and observed as S
        po[1,i,t,2] <- 0        # not yet entered and observed as I
        po[1,i,t,3] <- 1        # not yet entered and not observed
        po[2,i,t,1] <- pS[t]    # in S and observed as S
        po[2,i,t,2] <- 0        # in S and observed as I 
        po[2,i,t,3] <- 1-pS[t]  # in S and not observed 
        po[3,i,t,1] <- 0        # in I and observed as S
        po[3,i,t,2] <- pI[t]    # in I and observed as I
        po[3,i,t,3] <- 1-pI[t]  # in I and not observed
        po[4,i,t,1] <- 0        # dead and observed as S
        po[4,i,t,2] <- 0        # dead and observed as I
        po[4,i,t,3] <- 1        # dead and not observed
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
    
    
    # Calculate derived population parameters
    for (t in 1:(n.primary.occasions-1)){
      qgamma[t] <- 1-gammaS[t]-gammaI[t]  # prob of not entering
    }
    cprob[1] <- gammaS[1]+gammaI[1]   # cummulative prob of entering (either as S or I)
    for (t in 2:(n.primary.occasions-1)){
      cprob[t] <- (gammaS[t]+gammaI[t]) * prod(qgamma[1:(t-1)])
    } #t
    psi <- sum(cprob[])            # Inclusion probability
    for (t in 1:(n.primary.occasions-1)){
      b[t] <- (cprob[t] + 0.001) / (psi + 0.001)       # Entry probability
    } #t
    
    for (i in 1:nind){
      for (t in 1:n.primary.occasions){
        Ss[i,t] <- equals(z[i,t], 2) #
      } #t
      for (t in 1:n.primary.occasions){
        Is[i,t] <- equals(z[i,t], 3)
      } #t
      for (t in 2:(n.primary.occasions-1)){
        dS[i,t] <- equals(z[i,t]-Ss[i,t-1],0) 
        dI[i,t] <- equals(z[i,t]-Is[i,t-1],0)
      } #t   
      aliveS[i] <- sum(Ss[i,])
      aliveI[i] <- sum(Is[i,])
    } #i
    
    for (t in 1:n.primary.occasions){
      S[t]  <- sum(Ss[,t])       # Actual population size of S
      I[t]  <- sum(Is[,t])       # Actual pop size of I
      N[t]  <- S[t] + I[t]       # Actual total pop size
      BS[t] <- sum(dS[,t])       # Number of S entries
      BI[t] <- sum(dI[,t])       # Number of I entries
      B[t]  <- BS[t] + BI[t]     # total number of entries
      fS[t] <- (BS[t] + 0.001)/(N[t] + 0.001)# per capita recruitment rate of S
      fI[t] <- (BI[t] + 0.001)/(N[t] + 0.001)        # per capita recruitment rate of I
      f[t]  <- (B[t] + 0.001)/(N[t] + 0.001)         # total per capita recruitment rate
    } #t
    for (i in 1:nind){
      w[i] <- 1-equals(aliveS[i],0)-equals(aliveI[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
    
    }
    ",fill = TRUE)
sink()

