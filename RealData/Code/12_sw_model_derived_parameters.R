## Derived parameters ------------------------------------------------------- ##

 # 2/17/2020 - testing this code using the Zuni run from Zuni updated run
 # file is 19 GB - took ~855 to 859

 # make indexing easier - grabs everything up to the sims list so will have to
 # index what output (phi, alpha.0 etc.) plus index for that output (male, 
 # female, web etc.)
   model <- Z12.rCJS.MaxModel$BUGSoutput$sims.list

    
## MNA and Abundance Estimates ---------------------------------------------- ##


 # estimate MNA and abundance for each web


 # abundance - web 1
   for(m in 1:(n.primary.occasions)) {
     N[m] <- num.caught/(1 - (1 - p)^3)
   }


 # abundance - web 2
   for(m in 1:(n.primary.occasions)) {
     N[m] <- num.caught/(1 - (1 - p)^3)
   }


## -------------------------------------------------------------------------- ##  