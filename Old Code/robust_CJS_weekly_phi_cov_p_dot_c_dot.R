##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Modeling temporal and individual covariate, NDVI
#p and c are constant

#modeled on a weekly scale

#################################specify model in BUGS language
sink("robust_CJS_weekly_phi_ndvi_0_ndvi_1_tmax_3_tmin_5_p_dot_c_dot.bug")
cat("
    model{
    
    ###############Priors and constraints
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    alpha.ndvi_0 ~ dnorm(0, 0.4)T(-10,10)   #prior for slope on NDVI no lag
    alpha.ndvi_1 ~ dnorm(0, 0.4)T(-10,10)   #prior for slope on NDVI lag 1
    alpha.tmax_3 ~ dnorm(0, 0.4)T(-10,10)   #prior for slope on tmax lag 3
    alpha.tmin_5 ~ dnorm(0, 0.4)T(-10,10)   #prior for slope on tmin lag 5


    mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
    mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    

    for(i in 1:nind){
      for(w in 1:(n.weeks-1)){  
        # phi has  2 dimensions [indiv, and weeks]
        
        logit(phi[i,w]) <- alpha.0 + alpha.ndvi_0 * ndvi_0[i,w] + alpha.ndvi_1 * ndvi_1[i,w] + alpha.tmax_3 * tmax_3[i,w] + alpha.tmin_5 * tmin_5[i,w]

      } #w for weeks
    } #i for individual
    
    for(i in 1:nind){
      for(w in 1:n.weeks){  
        # p and c also have two dim 	
        logit(p[i, w]) <- mean.p  # could specify covariates here
        logit(c[i, w]) <- mean.c  
      } #w for weeks
    } #i for individual
    
    
    #############Likelihood 		
    # STATE PROCESS
    for(i in 1:nind){
      # define latent state at first capture 
      # dimensions [individual, week]
      z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
    
      for(w in (f[i]+1):n.weeks){  
        mu1[i, w] <- (phi[i, w-1]) * z[i, w-1] 
        z[i, w] ~ dbern(mu1[i, w]) 		
                
      } # w
    } # i
    
    # OBSERVATION PROCESS 
    for(obs in 1:n.obs){   
    
      y[obs] ~ dbern(z[id[obs], week[obs]] * ifelse(p.or.c[obs]==0,
        p[id[obs], week[obs]], c[id[obs], week[obs]]) ) 		
      # 0 represents p, 1 represents c (if caught before that session)



    } #obs
    
    } #model
    ",fill=TRUE)  
sink()
###########################################################################