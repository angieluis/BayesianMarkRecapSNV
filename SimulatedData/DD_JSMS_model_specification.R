## Jolly Seber Multistate Model specification Combined 3 sites using 
## Density Dependent K model. See SimulateDensityDependence_noinf_array.R 
## Jolly-Seber models allow for individuals to enter the population
## from a superpopulation - can estimate the prob of entrance
## See DD_JSMS_model_run.R for run code.



## To do: ----------------------------------------------------------##

## Error on indexing of N. Need to check all indexing.
## Check that the derived params and other indexing are accounting 
## for dummy occasion. (See book p325)
## And I don't want to calculate birth of prob of entrance from m1 to m2
## Because that assumes pop=0 at dummy occasion and would overestimate birth.
## Need priors for gamma[w,m]? Or change formulation? And need priors for that? 
## The book defines priors but doesn't give initial values or monitor them.
## Or how are they functions of derived birth rates?
## -----------------------------------------------------------------##

sink("DD_JSMS_model_specification.bug")
cat("
    model {

    # -------------------------------------------------
    # Parameters:
    # phi: survival probability 
    # p: recapture probability 
    # gamma: indiv prob of entry 
    # -------------------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    
    # Observations (O):
    # 1 seen 
    # 2 not seen
    # -------------------------------------------------

          
      ##### PRIORS FOR PHI #####
      m0          ~ dnorm(0, 1)T(-10, 10)    # prior for mortality intercept
      me          ~ dnorm(0, 0.4)T(-10, 10)    # prior for (logit) mortality at equil 
      m.male      ~ dnorm(0, 0.4)T(-10, 10)    # prior for additional male mortality

      #### PRIORS for birth ####
      b0          ~ dnorm(0, 0.4)T(-10, 10)   # prior for birth intercept (log at N=2)
      
      #### PRIORS for K  ######
      k.0      ~ dnorm(0, 1)T(-10, 10)    # prior for k coef on intercept for k
      k.ndvi      ~ dnorm(0, 0.4)T(-10, 10)    # prior for k coef on ndvi12
      

      ## seasonal covariates
      # don't estimate a beta coefficient for winter (assume 0),
      # then estimate a coefficent for all other seasons
      for(m in 1:3) {
        k.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on season
       }
      k.season.use <- c(0, k.season) 

	  
      ##### PRIORS FOR RECAPTURE #####
      sigma.0     ~ dnorm(0, 0.4)T(-10, 10)  # prior for intercept on p
      sigma.inf   ~ dnorm(0, 0.4)T(-10, 10)  # prior for infection on p





       
        ##### MODEL FOR PHI #####
      for(w in 1:n.webs){
        for(i in 1:n.inds[w]) {
          for(m in 1:(n.months[w] - 1)) {              
            
            # calculate K for that month based on ndvi and season
            log(K[w,m]) <- k.0  + k.ndvi*ndvi[w,m] + k.season.use[season[w,m]]
    
          #### Calculate web-specific N  
          
            Ns[i,m,w] <- equals(z[i,m,w], 2)  #  in book this had m-1 i think to remove dummy, but I changed

            N[w,m]  <- sum(Ns[,m,w])       # Actual population size  ---------------------------------------- how can it calc N when hasn't finished the m loop. Should i rearrange order of loops? or refer to m-1? or put outside of this set of loops? How do I get N to get z to get N? Circular. Indexing!

         
 
            ### Phi 
            logit(mort[i, m, w]) <-            
              m0 + (me-m0) * N[w,m]/K[w,m] +
              m.male * sex[i,w]
              
            phi[i, m, w] <- 1 - mort[i,m, w]

            } # m for months
          } # i for individual
      } # w for webs
        
        
       

        ##### MODEL FOR P #####
        ##### 3 dimensions [indiv, month, day, web]
     
      for(w in 1:n.webs){
        for(i in 1:n.inds[w]) {
          for(m in months.trapped.mat[w, 1:length.months.trapped[w]]) {

            # accounts for different secondary occasions
            for(d in 1:n.sec.occ) {               # simplified here for this model n.sec.occ doesn't vary, but could be n.sec.occ[w, Prim[i,m]
              logit(p[i, m, d, w]) <-
                sigma.0 


              } # d for days
            } # m for month
          } # i for individual
      } # w for web


 
  ######### Define state-transition and observation matrices --------------#
  for(w in 1:n.webs){
    for (i in 1:n.inds[w]){
      # Define probabilities of State(t+1) given State(t)
      for (m in 1:(n.months[w]-1)){
        ps[1,i,m,w,1] <- 1-gamma[m,w]               # still not yet entered
        ps[1,i,m,w,2] <- gamma[m,w]                 # just entered 
        ps[1,i,m,w,3] <- 0                          # not yet entered to dead
        ps[2,i,m,w,1] <- 0                          # alive to not yet entered
        ps[2,i,m,w,2] <- phi[i,m,w]                   # alive to alive (survival)
        ps[2,i,m,w,3] <- 1-phi[i,m,w]                # alive to dead
        ps[3,i,m,w,1] <- 0                          # dead to not yet entered
        ps[3,i,m,w,2] <- 0                          # dead to alive
        ps[3,i,m,w,3] <- 1                          # dead stay dead
      } #m
      
      for (m in months.trapped.mat[w, 1:length.months.trapped[w]]){        
        for(d in 1:n.sec.occ) {                           # ----------------------ditto simplified
        
          # Define probabilities of Observation(m) given State(m)
          # first index is state, last index is observation
          # could potentially include observation as wrong state (false neg or pos)
        po[1,i,m,d,w,1] <- 0              # not yet entered and observed 
        po[1,i,m,d,w,2] <- 1              # not yet entered and not observed
        po[2,i,m,d,w,1] <- p[i,m,d,w]     # alive and observed 
        po[2,i,m,d,w,2] <- 1-p[i,m,d,w]   # alive and not observed
        po[3,i,m,d,w,1] <- 0              # dead and observed 
        po[3,i,m,d,w,2] <- 1              # dead and not observed
        } # d
      } # m
     } # i
    } # w

    ############### Likelihood ---------------------------------------#
    # STATE PROCESS
    for(w in 1:n.webs){
      for (i in 1:n.inds[w]){
        # Define latent state at first dummy occasion
         z[i,1,w] <- 1   # Make sure that all individuals are in state 1 at t=1 (dummy occasion)
         # No one has entered yet (state 1) at t=1, because when input data above (Ch.primary.du), add a row before the actual data (where no has entered yet)

        for (m in 2:n.months[w]){  #
          # State process: draw S(t) given S(t-1)
          z[i,m,w] ~ dcat(ps[z[i,m-1,w], i, m-1, w, ])

        } #m
      } #i
    } #w

    # OBSERVATION PROCESS
    # draw O(t) given S(t)

  for(w in 1:n.webs){
    for (i in 1:n.inds[w]){
      # Define probabilities of State(t+1) given State(t)
      for (m in 2:(n.months[w])){
        for(d in 1:n.sec.occ){
        CHarray.aug.du[i, m, d, w] ~ dcat(po[z[i, m, w], i, m-1, d, w, ])  # think it should be m even tho book has m-1 
                 ## update this if not trapped every month (long.month, etc)
        } #d
      } #m
    } #i
  } #w





    ########################## Calculate derived population parameters
  for(w in 1:n.webs){  
    for (m in 1:(n.months[w]-1)){
      qgamma[m,w] <- 1-gamma[m,w] # prob of not entering
    }
    cprob[1,w] <- gamma[1,w]  # cummulative prob of entering 
    for (m in 2:(n.months[w]-1)){
      cprob[m,w] <- (gamma[m,w]) * prod(qgamma[1:(m-1),w])
    } #m
    psi <- sum(cprob[,w])            # Inclusion probability
    for (m in 1:(n.months[w]-1)){
      b[m,w] <- (cprob[m,w] + 0.001) / (psi + 0.001)      # Entry probability
    } #m
  } # w
  
    # for(w in 1:n.webs){             # don't need this because above for N in model
    #   for (i in 1:n.inds[w]){
    #     for (m in 2:n.months[w]){
    #       al[i,m-1,w] <- equals(z[i,m,w], 2)             
    #     } #t
    #     for (m in 1:(n.months[w]-1)){
    #       d[i,m,w] <- equals(z[i,m,w]-al[i,m,w],0)
    #     } #t   
    #     alive[i,w] <- sum(al[i, ,w])
    #   } #i
    # } # w
    
  be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
  
  for(w in 1:n.webs){
    for (m in 2:(n.months[w]-1)){
      B[m,w] <- sum(d[,m,w])       # Number entries   #the expected B[m,w]=M*cprob[m,w], where M is the augmented number of individuals (larger than pop)
      f[m,w]  <- (B[m,w] + 0.001)/(N[m,w] + 0.001)         # total per capita recruitment rate  
      
      f[m,w] ~ exp(b0 + (be-b0) * N[w,m-1]/K[w,m-1])  # think about the indexing here. Shoudl it be m-1? What about dummy occasion?
       ### is  ~ right here? 

    } #m
    for (i in 1:n.inds[w]){
      v[i,w] <- 1-equals(alive[i,w],0)
    } #i
    Nsuper[w] <- sum(v[,w])            # Superpopulation size per web
  }
}
    ",fill = TRUE)
sink()




