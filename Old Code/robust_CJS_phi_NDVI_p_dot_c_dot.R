##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Modeling temporal and individual covariate, NDVI
#p and c are constant

#################################specify model in BUGS language
sink("robust_CJS_phi_NDVI_p_dot_c_dot.bug")
cat("
    model{
    
    ###############Priors and constraints
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    alpha.1 ~ dnorm(0, 0.4)T(-10,10)   #prior for slope on NDVI
    mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
    mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    

    for(i in 1:nind){
      for(m in 1:(n.primary.occasions-1)){  
        # phi has  2 dimensions [indiv, and primary occasions]
        
        logit(phi[i,m]) <- alpha.0 + alpha.1 * NDVI[i,m] 

      } #m for months
    } #i for individual
    
    for(i in 1:nind){
      for(m in 1:n.primary.occasions){  
        # p and c also have two dim 	
        logit(p[i, m]) <- mean.p  # could specify covariates here
        logit(c[i, m]) <- mean.c  
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