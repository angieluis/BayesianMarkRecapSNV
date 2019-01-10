##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Modeling temporal and individual covariate, NDVI
#p and c are constant

#modeled on a weekly scale

#################################specify model in BUGS language
sink("robust_CJS_weekly_maxcov.bug")
cat("
    model{
    
    #########################Priors and constraints

    ### Priors for Phi
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    
    alpha.ndvi_1 ~ dnorm(0, 0.4)T(-10,10) #prior for coefficient on NDVI lag 1
    alpha.male ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on sex=male

    
      # don't estimate a beta coefficient for Jan (assume 0),
      # then estimate a coefficent for all other months 
    for(m in 1:11){
      alpha.month[m] ~ dnorm(0, 0.4)T(-10,10)   #prior for coefficients on months
      alpha.month.ndvi_1[m] ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on the interaction between month and ndvi lag 1
    }
    alpha.month.use <- c(0,alpha.month)
    alpha.month.ndvi_1.use <- c(0,alpha.month.ndvi_1)
    
    ###### Priors for recapture prob
    #mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
    #mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    sigma.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept on p
    sigma.recap ~ dnorm(0, 0.4)T(-10,10) #prior for coefficient on recapture in the the same occasion equivalent to c
    sigma.male ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on sex=male
    
    for(mo in 1:11){
      sigma.month[mo] ~ dnorm(0, 0.4)T(-10,10)   #prior for months
    }
    sigma.month.use <- c(0,sigma.month)


    # Model for Phi
    for(i in 1:nind){
      for(w in 1:(n.weeks-1)){  
        # phi has  2 dimensions [indiv, and weeks]
        

        logit(phi[i,w]) <- 
          alpha.0 + 
          alpha.ndvi_1 * ndvi_1[i,w] + 
          alpha.male * sex[i] +         #0 if female, 1 if male
          alpha.month.use[covariate.month[w]] +
          alpha.month.ndvi_1.use[covariate.month[w]] * ndvi_1[i,w]
      
      } #w for weeks
    } #i for individual
    
    # Model for p - 3 dimensions [indiv, week, day] 
    for(i in 1:nind){
      for(w in weeks.trapped){
        for(d in 1:n.sec.occ[covariate.prim[w]]){
          logit(p[i,w,d]) <- 
            sigma.0 +             # intercept
            sigma.recap * p.or.c[i,w,d] + #adjustment for if animal was caught previously in this primary session (0 if not caught before and 1 if so)
            sigma.male * sex[i] +   # adjustment for males (0 if female)
            sigma.month.use[covariate.month[w]] # month factor, where Jan=0
        } #d for days
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
    
      y[obs] ~ dbern(z[id[obs], week[obs]] * p[id[obs], week[obs], sec[obs]]) 		

    } #obs
    
    } #model
    ",fill=TRUE)  
sink()
###########################################################################