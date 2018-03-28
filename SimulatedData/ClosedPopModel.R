#######################################################################
#
# 6. Estimation of the size of a closed population
# 
#######################################################################

library(R2jags)

# 6.2. Generation and analysis of simulated data with data augmentation
# 6.2.1. Introduction to data augmentation for the simplest case: model M0
############################################### Simulate Data
# Define function to simulate data under M0
data.fn <- function(N = 100, p = 0.5, T = 3){
  yfull <- yobs <- array(NA, dim = c(N, T))
  for (j in 1:T){
    yfull[,j] <- rbinom(n = N, size = 1, prob = p)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs))
}


data <- data.fn()

str(data)

############################################ Analyze 

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
    model {
    
    # Priors
    omega ~ dunif(0, 1)   #prob of inclusion in the dataset (whether animal exists or not)
    p ~ dunif(0, 1)
    
    # Likelihood
    for (i in 1:M){
    z[i] ~ dbern(omega)			# Inclusion indicators
    for (j in 1:T){
    yaug[i,j] ~ dbern(p.eff[i,j])
    p.eff[i,j] <- z[i] * p		# Can only be detected if z=1
    } #j
    } #i
    
    # Derived quantities
    N <- sum(z[])
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 150))
abline(v = data$C, lwd = 3)





################################################ 6.2.2. Time effects: model Mt
# Define function to simulate data under Mt
data.fn <- function(N = 100, mean.p = 0.5, T = 3, time.eff = runif(T, -2, 2)){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.vec <- array(NA, dim = T)
  for (j in 1:T){
    p <- plogis(log(mean.p / (1-mean.p)) + time.eff[j])
    yfull[,j] <- rbinom(n = N, size = 1, prob = p)
    p.vec[j] <- p
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p.vec = p.vec, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn()

##############################
# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
    model {
    # Priors
    omega ~ dunif(0, 1)
    for (t in 1:T){
    p[t] ~ dunif(0, 1)
    }
    
    # Likelihood
    for (i in 1:M){   #i is individual
    z[i] ~ dbern(omega)
    for (t in 1:T){   #t is time
    yaug[i,t] ~ dbern(p.eff[i,t])
    p.eff[i,t] <- z[i] * p[t]   #mult by z so can't observe an animal that's not alive
    } #t
    } #i
    
    # Derived quantities
    N <- sum(z[])
    } # end model
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(data$T, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$BUGSoutput$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(70, 150))
abline(v = data$C, col = "black", lwd = 3)


########################### 6.2.3. Behavioral or memory effects: model Mb
# here c is when recaptured after caught in immediately previous session (if caught in 1 and not until 3, then both will be p's)

# Define function to simulate data under Mb
data.fn <- function(N = 200, T = 5, p = 0.3, c = 0.4){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.eff <- array(NA, dim = N)
  
  # First capture occasion
  yfull[,1] <- rbinom(n = N, size = 1, prob = p)
  
  # Later capture occasions
  for (j in 2:T){
    p.eff <- (1 - yfull[,(j-1)]) * p + yfull[,(j-1)] * c
    yfull[,j] <- rbinom(n = N, size = 1, prob = p.eff)
  }
  
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, c = c, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn(N = 200)

# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
    model {
    # Priors
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)     # Cap prob when caught at t-1
    c ~ dunif(0, 1)     # Cap prob when not caught at t-1
    
    # Likelihood
    for (i in 1:M){
    z[i] ~ dbern(omega)
    
    # First occasion
    yaug[i,1] ~ dbern(p.eff[i,1])
    p.eff[i,1] <- z[i] * p
    
    # All subsequent occasions
    for (j in 2:T){
    yaug[i,j] ~ dbern(p.eff[i,j])
    p.eff[i,j] <- z[i] * ( (1-yaug[i,(j-1)]) * p + yaug[i,(j-1)] * c )
    } #j
    } #i
    
    # Derived quantities
    N <- sum(z[])
    trap.response <- c - p
    } # end model
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "c", "trap.response", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)

hist(out$BUGSoutput$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(150, 300))
abline(v= data$C, col = "black", lwd = 3)




############################################ my modified version of behavioral model
# above formulation used a c only if caught in the immediately previous session
# we want c to be used if ever caught before (or caught before in the current session when we make robust design models). So p is first capture prob, and c is recapture prob. 

# 6.2.3. Behavioral or memory effects: model Mb
# Define function to simulate data under Mb
data.fn <- function(N = 200, T = 5, p = 0.3, c = 0.4){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.eff <- array(NA, dim = N)
  
  # First capture occasion
  yfull[,1] <- rbinom(n = N, size = 1, prob = p)
  
  # Later capture occasions
  for (j in 2:T){
    p.eff <- ifelse(sum(yaug[,1:(j-1)])==0,p,c) #if caught any time previously then c
    # will need to modify for robust design so if caught any time previously in this primary session
    yfull[,j] <- rbinom(n = N, size = 1, prob = p.eff)
  }
  
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, c = c, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn(N = 200)

# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model.jags")
cat("
    model {
    # Priors
    omega ~ dunif(0, 1)
    p ~ dunif(0, 1)     # Cap prob when never caught before
    c ~ dunif(0, 1)     # Cap prob when caught before
    
    # Likelihood
    for (i in 1:M){
    z[i] ~ dbern(omega)
    
    # First occasion
    yaug[i,1] ~ dbern(p.eff[i,1])
    p.eff[i,1] <- z[i] * p
    
    # All subsequent occasions
    for (t in 2:T){
    yaug[i,t] ~ dbern(p.eff[i,t])
    p.eff[i,t] <- z[i] * ifelse(sum(yaug[i,1:(t-1)])==0,p,c)
    } #t
    } #i
    
    # Derived quantities
    N <- sum(z[])
    trap.response <- c - p
    } # end model
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "c", "trap.response", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT <1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)

hist(out$BUGSoutput$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(150, 300))
abline(v= data$C, col = "black", lwd = 3)

### quite often  overestimates