##### Robust Design CJS with capture histories as array [indiv,prim,sec]

# Modeling temporal and individual covariate, NDVI
#p and c are constant

#modeled on a monthly scale

## in the process of updating code to account for months not trapped at specific webs - update p code here particular.

#################################specify model in BUGS language
sink("robust_CJS_monthlylist_CovInd.bug")
cat("
    model{
    
    #########################Priors and constraints

    ### Priors for Phi
    alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
    
    alpha.male ~ dnorm(0, 0.4)T(-10,10)   #prior for coef on sex=male

    
      # don't estimate a beta coefficient for Jan (assume 0),
      # then estimate a coefficent for all other months 
    for(m in 1:11){
      alpha.month[m] ~ dnorm(0, 0.4)T(-10,10)   #prior for coefficients on months
      
    }
    alpha.month.use <- c(0,alpha.month)
    
    ## add indicator variable for 3, 6 , 12 month lag for NDVI

    for(i in 1:n.lags){
      ind[i] ~ dbern(0.2)
      lag.coefT[i] ~ dnorm(0,0.4)
      lag.coef[i] <- ind[i]*lag.coefT[i]
    }




    
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
      for(m in 1:(n.months-1)){  
        # phi has  2 dimensions [indiv, and months]
        

        logit(phi[i,m]) <- 
          alpha.0 + 
          alpha.male * sex[i] +         #0 if female, 1 if male
          alpha.month.use[covariate.month[m]] +
          inprod(lag.mat[i,],lag.coef)
      } #m for months
    } #i for individual
    
    # Model for p - 3 dimensions [indiv, month, day] #same number of months as phi (longmonths), but with NAs for months/days not trapped
    for(i in 1:nind){
      for(m in months.trapped.mat[i, 1:length.months.trapped[i]]){
        for(d in 1:n.sec.occ[Prim[m]]){
          logit(p[i,m,d]) <- 
            sigma.0 +             # intercept
            sigma.recap * p.or.c[i,m,d] + #adjustment for if animal was caught previously in this primary session (0 if not caught before and 1 if so)
            sigma.male * sex[i] +   # adjustment for males (0 if female)
            sigma.month.use[covariate.month[m]] # month factor, where Jan=0
          
        } #d for days
      } #m for month
    } #i for individual
    
    
    #############Likelihood 		
    # STATE PROCESS
    for(i in 1:nind){
      # define latent state at first capture 
      # dimensions [individual, month] spans study period not just months trapped
      z[i,f.longmonth[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
    
      for(m in (f.longmonth[i]+1):n.months){  
        mu1[i, m] <- phi[i, m-1] * z[i, m-1] 
        z[i, m] ~ dbern(mu1[i, m]) 		
                
      } # m (total months spanned)
    } #i
    
    # OBSERVATION PROCESS 
    for(obs in 1:n.obs){   
      
      y[obs] ~ dbern(z[id[obs], longmonth.obs[obs]] * p[id[obs], longmonth.obs[obs], sec[obs]]) 		
      
      
    } #obs
    

    } #model
    ",fill=TRUE)  
sink()
###########################################################################