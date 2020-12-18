##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Model assuming phi is monthly taking into account the monthly covariates for months not trapped 
#p and c are constant

# This ignores p and c for first capture occasion. So both phi and p vectors are of length of primary sessions minus one. Gets put back in if using Jolly Seber and estimating N.

#################################specify model in BUGS language
sink("robust_CJS_phi_month_p_dot_c_dot2.bug")
cat("
    model{
    
    ###############Priors and constraints
    for(u in 1:12){
      mean.phi[u] ~ dnorm(0, 0.4)T(-10,10)     # priors for monthly survival
    }
    mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
    mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    
    for(i in 1:nind){
      for(m in 1:(n.primary.occasions-1)){  
        # phi has  2 dimensions [indiv, and primary occasions]
        
        months <- long.month[which(covariate.prim==m)] 
        # if there are no gaps, months should be length 1, but if some months 
        # weren't trapped there may be multiple months' data which should be 
        # incorporated
##### BUGS isn't recognizing the function which - can I say long.month[covariate.prim==m]

        denom <- 1 + exp(-mean.phi[covariate.month[which(covariate.prim==m)]])
    
        phi[i,m] <- (1/prod(denom))^(1/length(months))
    
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