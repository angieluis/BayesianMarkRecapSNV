logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}


# function to create a primary CH from the secondary capture history:
primaryMS.fun<-function(CH.secondary){
  x <- lapply(CH.secondary,function(x){apply(x,1,min)})
  v1 <- unlist(x)
  CH.primary <- matrix(v1, nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary))
  return(CH.primary)
}

#functions to add dummy occasion 
primary.dummy.fun <- function(CH.primary,notseen=3){
  CH.primary.du <- cbind(rep(notseen, dim(CH.primary)[1]), CH.primary)
  return(CH.primary.du)
}

secondary.dummy.fun <- function(CH.secondary,notseen=3){
  CH.secondary.du <- c(list(matrix(notseen,nrow=dim(CH.secondary[[1]])[1],ncol=dim(CH.secondary[[1]])[2])), CH.secondary)
  return(CH.secondary.du)
}

# functions to Augment data
primary.augment.fun <- function(CH.primary.du,notseen=3,num.aug=500){
  nz <- num.aug
  CH.primary.ms <- rbind(CH.primary.du, matrix(notseen, ncol = dim(CH.primary.du)[2], nrow = nz))
  return(CH.primary.ms)
}

secondary.augment.fun <- function(CH.secondary.du,notseen=3,num.aug=500){
  nz <- num.aug
  CH.secondary.ms <- lapply(CH.secondary.du,function(x){rbind(x,matrix(notseen, ncol = dim(x)[2], nrow = nz))})
  return(CH.secondary.ms)
}


# Function to create known latent states z
### fill in all known but unobserved states (can't go back to S from I)
#If observed as I, then not seen, and seen again later, when not seen must have been I.
#If observed as S, not seen, then observed as S again, then must be S.
# Remember these are now states not observations, so coded differently.
# 1 is not yet entered
# 2 is S
# 3 is I
# 4 is dead
# only able to fill in 2's and 3's
# Allows us to fill in a lot. And should speed up computation time
known.state.SImsJS <- function(ms=CH.primary.ms, notseen=3){ # ms is multistate capture history
  # notseen: label for 'not seen' #here is 3
  state <- ms
  state[state==notseen] <- NA
  for(i in 1:dim(ms)[1]){
    if(length(which(ms[i, ] == 2)) > 0){ #filling in I's where can
      minI <- min(which(ms[i, ] == 2)) #I's are observation 2
      maxI <- max(which(ms[i, ] == 2))
      state[i, minI:maxI] <- 3}         # I's are state 3
    if(length(which(ms[i, ]==1)) > 0){  #filling in S's where can
      minS <- min(which(ms[i, ] == 1))  # S's are observation 1
      maxS <- max(which(ms[i, ] == 1))
      state[i, minS:maxS] <- 2}         # S's are state 2
  }
  return(state)
}



# Specify initial values
jsmsinf.init <- function(ch=CH.primary.ms, num.aug=500){
  # ch is primary capture histories after augmentation
  # nz is number of rows added for augmentation
  nz <- num.aug
  kn.state <- known.state.SImsJS(ms=ch)
  state <- matrix(2, nrow=dim(ch)[1], ncol=dim(ch)[2]) # default is S (2)
  state <- replace(state,!is.na(kn.state),NA)
  
  for(i in 1:(dim(state)[1]-nz)){
    f <- min(which(is.na(state[i,])))       # before ever caught
    if(f>1){state[i,1:(f-1)] <- 2}              # tried both 1 and 2 here, still get errors
    if(length(which(kn.state[i,] == 3)) > 0){
      maxI <- max(which(kn.state[i,]==3))
      if(maxI<dim(state)[2] ){
        state[i,(maxI+1):dim(state)[2]] <- 3 # all after caught as I are I (3)
      }
    }
  }
  state[(dim(state)[1]-nz+1):dim(state)[1],] <- 1
  state[,1] <- NA #this is specified in likelihood
  return(state)
}

