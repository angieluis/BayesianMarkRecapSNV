##################################################################################
## Robust design Multistate: seronegative, seropositive, (unknown serostatus)
##  measuring survival (phi), recapture (p), probabilty of becoming infected (Psi)
##################################################################################

#-----------------To Do------------------------------------------#
# 1) Deal with unknown serostatus: put in NA's in y_infected
#     Still need to figure out the last few lines of code
# 2) transforms/priors for psi/beta
# 3) check the indexing
# 4) make this line up with newest code for inconsistent trapping
# 5) Make sure all priors are listed
#  (Working with "Z12cut_RobustMS.R" model run code)
#----------------------------------------------------------------#


# Specify model in BUGS language
sink("RDMS_phi_sexinf2.bug")
#cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability with dimensions [i,m,inf]
    # psiSI:  probability of becoming infected (from S to I)
    #     now a function of beta*I (proper process model) 
    # p: recapture probability with dimensions [i,m,d,inf]
      
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
    
      
        psiSI[i,t] <- beta[i]*I.dat[i,t] # /N.dat[i,t]^q  # think about indexing here 
    # really should calculate I here based on p (num.caught/p). but I think there will be estimation problems, so for now, I'll go with using MNI (scaled between (0,1))
    
    # I don't want a logit transform on beta/psi because of the S shape
    # so need to make sure beta is bounded so that beta*max(I.dat)<1.
    # Since I.dat will be scaled with max 1, beta needs to be between 0 and 1.
    # This also doens't account for infected immigration. That means if I[t-1]=0,
    # then no chance of getting infected.  
    # Could also change to mean.beta*I[t-1]/N^q and estimate q to account for 
    # frequency dependent transmission or in between. But then need to worry about the bounds 
    
      } #i
    } #t
    
    
    for (u in 1:2){     #for both states
      sigma.0[u] ~ dnorm(0, 0.4)T(-10,10)      # Priors for mean state-spec. recapture
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
    
    for (t in 1:(n.primary.occasions-1)){
      for(i in 1:nind){
        for(inf in 1:2){
          logit(phi[i,t,inf]) <- 
            alpha.0[inf] +                       # intercept
            alpha.male * sex[i]           # sex
          # might want to add interaction between infection and sex
        } # inf
        
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of State(t+1) given State(t)
      for (t in 1:(n.primary.occasions-1)){
        ps[1,i,t,1] <- phi[i,t,1] * (1-psiSI[i,t]) # survival of S to S
        ps[1,i,t,2] <- phi[i,t,1] * psiSI[i,t]     # survival and transition from S to I
        ps[1,i,t,3] <- 1-phi[i,t,1]              # S to dead
        ps[2,i,t,1] <- 0                      # I to S (can't happen, tho check data)
        ps[2,i,t,2] <- phi[i,t,2]                # survival of I to I
        ps[2,i,t,3] <- 1-phi[i,t,2]              # I to dead
        ps[3,i,t,1] <- 0                      # dead to S
        ps[3,i,t,2] <- 0                      # dead to I
        ps[3,i,t,3] <- 1                      # dead stay dead
      }
      
      
      ### make p a 4 dimensional array that is p[i,m,d,inf]
      # where [i,m,d] refer to individual, temporal covariates and p.or.c
      # and inf dimension has length 2, where 1=uninfected, and 2=infected
      # so doesn't strictly refer to individual p, because an infected individual will
      # still have an entry for uninfected dimension
      for(i in 1:nind){
        for(m in months.trapped.mat[i, 1:length.months.trapped[i]]){
          for(d in 1:n.sec.occ[Prim[m]]){
            for(inf in 1:2){
              logit(p[i,m,d,inf]) <- 
                sigma.0[inf] +             # intercept
                sigma.recap * p.or.c[i,m,d] + #adjustment for if animal was caught previously in this primary session (0 if not caught before and 1 if so)
                sigma.male * sex[i] +         # adjustment for males (0 if female)
                sigma.season.use[season[m]] + # month factor, where Jan=0
                sigma.year.use[year[m]] +     # year factor, where first year=0
                sigma.web.use[web[i]]
            } # inf for infection 1 or 2
          } #d for days
        } #m for month
      
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
        
      y_caught[obs] ~ dbern(z[id[obs], longmonth.obs[obs]] * p[id[obs], longmonth.obs[obs], sec[obs], ifelse(z[id[obs], longmonth.obs[obs]]==2,2,1)])
      
      y_infected[obs] ~ y_caught[obs] * ifelse(z[id[obs], longmonth.obs[obs]]==2,1,0)
        # this isn't right because not a random draw. should i just change to equal?
        # don't think I can.
    } #obs
    
    
    
    }
    ",fill = TRUE)
sink()

