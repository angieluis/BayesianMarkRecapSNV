##### Robust Design CJS with capture histories as list and manipulated (see RobustCJSRaggedArrayParallel.R)

# Modeling temporal and individual covariate, NDVI
#p and c are constant

#modeled on a weekly scale

#################################specify model in BUGS language
sink("robust_CJS_weekly_phi_ndvi_0_ndvi_1_tmax_3_tmin_5_p_dot_c_dot.bug")
cat("
    model{
    
    #########################Priors and constraints

    ### Priors for Phi
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    
    alpha.ndvi_1 ~ dnorm(0, 0.4)T(-10,10) #prior for coefficient on NDVI lag 1
    alpha.male ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on sex=male
    
    # for month covariate, what are all the months in the dataset 
    # (e.g, may not ever have trapped in Jan)
    uniq.months <- sort(unique(covariate.month))
      # don't estimate a beta coefficient for the first month (assume 0),
      # then estimate a coefficent for all other months (length(uniq.months)-1)
    for(m in 1:(length(uniq.months)-1)){
      alpha.month[m] ~ dnorm(0, 0.4)T(-10,10)   #prior for coefficients on months
      alpha.month.ndvi_1[m] ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on the interaction between month and ndvi lag 1
    }
    
    ###### Priors for recapture prob
    #mean.p ~ dnorm(0, 0.4)T(-10,10)       # prior for p
    #mean.c ~ dnorm(0, 0.4)T(-10,10)       # prior for c
    sigma.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    
    sigma.recap ~ dnorm(0, 0.4)T(-10,10) #prior for coefficient on recapture in the the same occasion (trap happy/trap shy) equivalent to c
    sigma.male ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on sex=male
    for(m in 1:(length(uniq.months)-1)){
      sigma.month[m] ~ dnorm(0, 0.4)T(-10,10)   #priors for months
    }
    
    # Model for Phi
    for(i in 1:nind){
      for(w in 1:(n.weeks-1)){  
        # phi has  2 dimensions [indiv, and weeks]
        cmonth <- covariate.month[w] # current month

        logit(phi[i,w]) <- 
          alpha.0 + 
          alpha.ndvi_0 * ndvi_0[i,w] + 
          alpha.ndvi_1 * ndvi_1[i,w] + 
          alpha.male * sex[i] +         #0 if female, 1 if male
          ifelse(cmonth==uniq.months[1],0, alpha.month[which(uniq.months==cmonth)-1]) +
          ifelse(cmonth==uniq.months[1],0, alpha.month.ndvi_1[which(uniq.months==cmonth)-1]) * ndvi_1[i,w]

      } #w for weeks
    } #i for individual
    
    # Model for p
    for(i in 1:nind){
      for(w in 1:n.weeks){ 
        for(d in 1:n.sec.occ[w]){
        # p has 3 dimensions [indiv, week, day] 	
          logit(p[i,w,d]) <- 
            sigma.0 +             # intercept
            sigma.recap * p.or.c[i,w,d] + #adjustment for if animal was caught previously in this primary session (0 if not caught before and 1 if so)
            sigma.male + sex[i] +   # adjustment for males (0 if female)
            ifelse(cmonth==uniq.months[1],0, sigma.month[which(uniq.months==cmonth)-1]) # month factor, where first unique month=0
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