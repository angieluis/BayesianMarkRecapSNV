## Model - Robust design CJS Covariate Indicator ---------------------------- ##


 # capture histories as array [indiv, prim, sec]
 # indicator estimates probability that a cov. should be included in the model
 # https://darrenjw.wordpress.com/2012/11/20/getting-started-with-bayesian-variable-selection-using-jags-and-rjags/


## Specify model ------------------------------------------------------------ ##


 # divert - can rename based on model
 # updated 11/19/2019 - removed 18-month time lag from 
   sink("Code/CovInd_Interaction_Updated.bug")
   
 
 # print to file
   cat(
   
     model {
       
       ##### PRIORS AND CONSTRAINTS #####
       
       ##### PRIORS FOR PHI
       alpha.0 ~ dnorm(0, 0.4)T(-10,10)   #prior for intercept
       alpha.male ~ dnorm(0, 0.4)T(-10, 10)   # prior for coeff on sex = male
       
       
       # don't estimate a beta coefficient for winter (assume 0 = winter, then 
       # spring, summer, fall) then estimate a coefficent for all other seasons
       for(m in 1:3) {
         alpha.season[m] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coeff on months
         }
       alpha.season.use <- c(0, alpha.season)
       
       
       # coefficient for webs, 0 for first web etc.
       for(w in 1:(max(web)-1)) {
         alpha.web[w] ~ dnorm(0, 0.4)T(-10, 10)   # prior for coeff on web
         }
       alpha.web.use <- c(0, alpha.web)
       
       
       # beta distribution suggested - see below link
       # pind ~ dbeta(1, 3) # prior on proportion of covariates to include
       
       
       # add indicator variable for variables in the array
       for(cov in 1:n.covariates) {
         # use p.ind (probability that we should include the covariate)
         # bernoulli distribution draws a 1 with that probability otherwise a 0
         # ind[cov] ~ dbern(pind)
         ind[cov] ~ dbern(1/n.covariates)
            # make a random draw (0 or 1)
            # 1/number of covariates
         for(m in 1:4){
            cov.coefT[cov,m] ~ dnorm(0, 0.4)  
            cov.coef[cov,m] <- ind[cov] * cov.coefT[cov, m]
         }
       }
       
       
       ##### PRIORS FOR RECAPTURE
       sigma.0 ~ dnorm(0, 0.4)T(-10, 10)   # prior for intercept on p
    
       # multiplies sigma recap times the p.or.c.array
       # prior for coefficient on recapture in the same occasion, equiv. to c
       sigma.recap ~ dnorm(0, 0.4)T(-10, 10) 
       sigma.male ~ dnorm(0, 0.4)T(-10, 10)   # prior for coef on sex=male
       
       
       for(se in 1:3) {
         sigma.season[se] ~ dnorm(0, 0.4)T(-10, 10)   #prior for months
         }
       sigma.season.use <- c(0,sigma.season)
       
       
       # prior for coeff on web
       for(w in 1:(max(web)-1)) {
         sigma.web[w] ~ dnorm(0, 0.4)T(-10, 10)   #prior for coef on web
         }
       sigma.web.use <- c(0, sigma.web)
       
       
       ##### MODEL FOR PHI #####
       for(i in 1:nind) {
         for(m in 1:(n.months - 1)) {
           
           # phi has  2 dimensions [indiv, and months]
           logit(phi[i,m]) <- 
             alpha.0 + 
             alpha.male * sex[i] +         #0 if female, 1 if male
             alpha.season.use[season[m]] +  
             alpha.web.use[web[i]] +
               
              
             # for interaction between season and covariate
             inprod(covariate.array[i, m, ], cov.coef[ , season[m]])
           } # m for months
         } # i for individual
       
       
       ##### MODEL FOR P #####
       # 3 dimensions [indiv, month, day]
       # same number of months as phi (longmonths), NA for times not trapped
       for(i in 1:nind) {
         for(m in months.trapped.mat[i, 1:length.months.trapped[i]]) {
           for(d in 1:n.sec.occ[i,Prim[m]]) {
             logit(p[i,m,d]) <- 
               sigma.0 +             # intercept
               sigma.recap * p.or.c[i,m,d] +   # adjustment for if animal was  
                                               # caught previously in this 
                                               # primary session (0 if not 
                                               # caught before and 1 if so)
               sigma.male * sex[i] +           # male adjustment (0 if female)
               sigma.season.use[season[m]] +   # season factor, where winter=0
               sigma.web.use[web[i]]
             } # d for days
           } # m for month
         } # i for individual
       
       
       ##### LIKELIHOOD #####
       
       ##### STATE PROCESS
       for(i in 1:nind) {
         # define latent state at first capture 
         # dimensions [individual, month] spans study period, all months
         # z is true (latent) state alive or dead, know alive at first capture
         z[i,f.longmonth[i]] <- 1
         
         for(m in (f.longmonth[i] + 1):n.months) {  
           mu1[i, m] <- phi[i, m - 1] * z[i, m - 1]
           z[i, m] ~ dbern(mu1[i, m])
           } # m (total months spanned)
         } #i
       
       ##### OBSERVATION PROCESS 
       for(obs in 1:n.obs) {
         y[obs] ~ dbern(z[id[obs], 
                          longmonth.obs[obs]] * p[id[obs], 
                                                  longmonth.obs[obs], 
                                                  sec[obs]])
         } # obs
       } # model
     , fill = TRUE)


 # sink file
   sink()


## -------------------------------------------------------------------------- ##