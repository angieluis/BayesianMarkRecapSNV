##### Functions

logit=function(x){
  log(x/(1-x))}
revlogit=function(x){
  exp(x)/(1+exp(x))}

# function to create primary CH from secondary list
primary.ch.fun <- function(CH.secondary){ # as list of monthly matrices 
  
  x <- lapply(CH.secondary,rowSums)
  v1 <- unlist(x)
  CH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 
  
  return(CH.primary)
}


###### function to create matrix of initial values for latent state z
# we shouldn't give initial values for those elements of z whose value is specified in the data.
# they get an NA
cjs.init.z=function(ch,f){
  for(i in 1:dim(ch)[1]){
    if(sum(ch[i,])==1) next
    n2=max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)	
}

#function to create matrix with info about known latent state z
known.state.cjs=function(ch){
  state=ch
  for(i in 1:dim(ch)[1]){
    n1=min(which(ch[i,]==1))
    n2=max(which(ch[i,]==1))
    state[i,n1:n2]=1
    state[i,n1]=NA			#only filling in those that were 0s but we know were alive because caught before and after
  }
  state[state==0]=NA
  return(state)
}





#### add months not trapped so can have temporal covariates for all months, and line them up to sessions. so there will be multiple temporal covariates assigned to one primary session if months not trapped.

# This applies to survival (phi) 
# if want to apply covariates to capture probabilities, then months not trapped won't matter, only use the first one

temporaldata.fun <-function(data,site,web=NULL){
  if(length(web)==1){
    sessions <- sort(unique(data$Session[which(data$site==site&data$web==web)]))
    }else{
      sessions <- sort(unique(data$Session[which(data$site==site)]))
    }
  
  first.session <- sessions[1]
  last.session <- sessions[length(sessions)]
  first.montha <- strsplit(as.character(first.session),split=character(0))[[1]][5:6]
  first.month <- as.numeric(paste(first.montha[1],first.montha[2],sep=""))
  
  first.yeara <- strsplit(as.character(first.session),split=character(0))[[1]][1:4]
  first.year <- as.numeric(paste(first.yeara[1],first.yeara[2],first.yeara[3],first.yeara[4],sep=""))
  last.montha <- strsplit(as.character(last.session),split=character(0))[[1]][5:6]
  last.month <- as.numeric(paste(last.montha[1],last.montha[2],sep=""))
  last.yeara <- strsplit(as.character(last.session),split=character(0))[[1]][1:4]
  last.year <- as.numeric(paste(last.yeara[1],last.yeara[2],last.yeara[3],last.yeara[4],sep=""))
  
  y <- last.year-first.year
  ms <- c(first.month:12,rep(1:12,y-1),1:last.month)
  ys <- c(rep(first.year,length(first.month:12)),rep((first.year+1):(last.year-1),each=12),rep(last.year,length(1:last.month)))
  yearmonth <- paste(as.character(ms),as.character(ys),sep="")
  mc <- as.character(ms)
  for(i in 1:length(ms)){
    mc[i] <- ifelse(ms[i]<10,paste(paste(as.character(ys[i]),"0",as.character(ms[i]),sep="")),paste(paste(as.character(ys[i]),as.character(ms[i]),sep="")))
  }
  
  s1 <- sessions[match(mc,sessions)]
  not.trapped <- which(is.na(s1))
  #mc[not.trapped]
  
  sn <- 1:length(sessions)
  session.num <- sn[match(mc,sessions)]
  for(i in 1:length(not.trapped)){
    session.num[not.trapped[i]]<-session.num[max(which(sn<not.trapped[i]))]
  }
  
  temporal.covariates <- data.frame(long.month=1:length(ms),session= s1,year=ys,month=ms,covariate.prim=session.num)
  
  return(temporal.covariates)
}



###############################
## Create a dataframe for individual covariates

####### Individual Covariates

individual.covariate.fun <- function(data, tags, Ch.secondary){
  ic <- data.frame(ID=1:dim(Ch.secondary[[1]])[1],tag=tags)

  web <- character()
  sex <- numeric() #1 male, 0 female
  for(i in 1:dim(Ch.secondary[[1]])[1]){
    ind <- which(data$tag==IDs[i])
    x <- data[ind,]
    x <- x[order(x$Session),]
    web[i] <- as.character(x$web[1])
    
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- web
  ic$sex <- replace(sex,sex==-9,NA)
  return(ic)
}