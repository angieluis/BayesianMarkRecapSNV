
  model {
  #
  # RD model phi(.) p(.)
  #
  # priors
    b0.phi ~ dnorm( 0, 0.001 )T(-10,10)
    b0.p ~ dnorm( 0, 0.001 )T(-10,10)
  
    logit(mean.phi) <- b0.phi
    logit(mean.p) <- b0.p
  
    for( i in 1:nind ){
      for( t in f[i]:(nocc-1) ){
          phi[i,t] <- mean.phi
      }
    }
    for( i in 1:nind ){
      for( t in f[i]:(nocc) ){
        for( tSec in 1:nSecOcc[t] ){
          p[i,t,tSec] <- mean.p
        }
      }
    }
  
    # likelihood
    for( i in 1:nind ){
      z[i,f[i]] <- 1
      for( tSec in 1:nSecOcc[f[i]] ){
        y[i,f[i],tSec] ~ dbern( mu2[i,f[i],tSec] )
        mu2[i,f[i],tSec] <- p[i,f[i],tSec] * z[i,f[i]]
      }
      for( t in (f[i]+1):nocc ){
      
        # state
        z[i,t] ~ dbern( mu1[i,t] )
        mu1[i,t] <- phi[i, t-1] * z[i,t-1]
        
        # observation
        for( tSec in 1:nSecOcc[t] ){
          y[i,t,tSec] ~ dbern( mu2[i,t,tSec] )
          mu2[i,t,tSec] <- p[i,t,tSec] * z[i,t]
        }
      }
    }
  
  }
      