##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Model assuming phi, p and c are a function of month (factors so 1 estimate per month for each)

# This ignores p and c for first capture occasion. So both phi and p vectors are of length of primary sessions minus one. Gets put back in if using Jolly Seber and estimating N.

#################################specify model in BUGS language
sink("robust_CJS_phi_month_p_month_c_month.bug")
cat("
    model{
    
    ###############Priors and constraints
    for(u in 1:12){
      mean.phi[u] ~ dnorm(0, 0.4)T(-10,10)     # priors for monthly survival
      mean.p[u] ~ dnorm(0, 0.4)T(-10,10)       # prior for p
      mean.c[u] ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    }

    for(i in 1:nind){
      for(m in 1:(n.primary.occasions-1)){  
        # phi has only 2 dimensions [indiv, and primary occasions]
        logit(phi[i,m]) <- mean.phi[month[m]]   # could specify covariates here
            # this is weekly survival. 
            # Below converted to survival over whole primary period (usually 4-5 weeks).
            # For time-varying covariate when months weren't sampled need to put in 
            # conditional statement here so that I include all the covariates 
            # for all the months included
      } #m for months
    } #i for individual
    
    for(i in 1:nind){
      for(m in 1:n.primary.occasions){  
        # p and c also have two dim (assuming won't vary by day within session, e.g. weather)	
        logit(p[i, m]) <- mean.p[month[m]]  # could specify covariates here
        logit(c[i, m]) <- mean.c[month[m]]  
      } #m for months
    } #i for individual
    
    
    #############Likelihood 		
    # STATE PROCESS
    for(i in 1:nind){
      # define latent state at first capture 
      # dimensions [individual, primary session (month)]
      z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
    
      for(m in (f[i]+1):n.primary.occasions){  
        mu1[i, m] <- (phi[i, m-1])^time.int[m-1] * z[i, m-1] 
        z[i, m] ~ dbern(mu1[i, m]) 		
                
      } # m
    } # i
    
    # OBSERVATION PROCESS 
    for(obs in 1:n.obs){   
    
      y[obs] ~ dbern(z[id[obs], prim[obs]] * ifelse(p.or.c[obs]==0, p[id[obs], prim[obs]],            c[id[obs], prim[obs]]) ) 		# 0 represents p, 1 represents c (if caught before that       session)

    } #obs
    
    } #model
    ",fill=TRUE)  
sink()
###########################################################################