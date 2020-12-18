### Functions to combine data from multiple sites (all these will work for multistate models as well as CJS)
#############################################################

library(rlist)




# paste all the observation data together (in alphabetical order)
# change IDs so match with individual covariates above.
combine.obsdat.fun <- function(obs.dats = list(obs.dat_GC,obs.dat_N,obs.dat_Z), 
            individual.covariates = list(covariate.data_GC$individual.covariates,covariate.data_N$individual.covariates, covariate.data_Z$individual.covariates)){

  n.sites <- length(individual.covariates)
  site.names <- unlist(lapply(individual.covariates,function(x){unlist(strsplit(as.character(x$tag[1]),"[.]"))[1]}))
  num.IDs <- unlist(lapply(individual.covariates,function(x){dim(x)[1]})) # number of individuals in data sets
  
  IDs <- list(data.frame(new=individual.covariates[[1]]$ID,old=individual.covariates[[1]]$ID))
  for(i in 2:n.sites){
    old <- individual.covariates[[i]]$ID
    s <- IDs[[i-1]]$new[length(IDs[[i-1]]$new)]+1
    new <- s:(s+num.IDs[i]-1)
    IDs[[i]] <- data.frame(new,old)
  }
  
  #renumber obs.dat with sequential numbering
  for(i in 2:n.sites){
    obs.dats[[i]]$ID <- IDs[[i]]$new[match(obs.dats[[i]]$ID,IDs[[i]]$old)]
  }
  
  
  combined.obsdat <- obs.dats[[1]]
  for(i in 2:n.sites){
    combined.obsdat <- rbind(combined.obsdat,obs.dats[[i]])
  }
  return(combined.obsdat)
}





# combine all covariates 
combine.covariates.fun <- function(covariate.data = list(covariate.data_GC,covariate.data_N, covariate.data_Z)){
  
  n.sites <- length(covariate.data)
  site.names <- unlist(lapply(covariate.data ,function(x){unlist(strsplit(as.character(x$individual.covariates$tag[1]),"[.]"))[1]}))
  num.IDs <- unlist(lapply(covariate.data,function(x){dim(x$individual.covariates)[1]})) # number of individuals in data sets
  times <- unlist(lapply(covariate.data,function(x){dim(x$temporal.covariates)[1]}))
  max.time <- max(times)
  num.covariates <- length(covariate.data[[1]])-2 ## assumes that the first two elements are the individual and temporal covariate data frames and the rest are covariate matrices
  
  IDsites <- data.frame(site=site.names,start=c(1,rep(NA,n.sites-1)), end=c(num.IDs[1],rep(NA,n.sites-1)))
  for(i in 2:n.sites){
    IDsites$start[i] <- IDsites$end[i-1]+1
    IDsites$end[i] <- IDsites$start[i] + num.IDs[i]-1
  }
  
  # add site and max month covariate
  for(i in 1:n.sites){
    covariate.data[[i]]$individual.covariates$site <- site.names[i]
    covariate.data[[i]]$individual.covariates$max.months <- times[i]
  }
  
  #change ID so sequential
  for(i in 2:n.sites){
    s <- covariate.data[[i-1]]$individual.covariates$ID[length(covariate.data[[i-1]]$individual.covariates$ID)]+1
    covariate.data[[i]]$individual.covariates$ID <- s:(s+num.IDs[i]-1)
  }
  
    
  combined.individual.cov.data <- covariate.data[[1]]$individual.covariates
  for(i in 2:n.sites){
    combined.individual.cov.data <- rbind(combined.individual.cov.data,covariate.data[[i]]$individual.covariates)
  }
  
  covariate.matrices <- list()
  #for each covariate paste together 
  for(i in 1:num.covariates){
    
    covariate.matrices[[i]] <- matrix(NA, nrow=sum(num.IDs), ncol=max.time)
    
    for(j in 1:n.sites){
    
      covariate.matrices[[i]][IDsites$start[j]:IDsites$end[j], 1:times[j]] <- covariate.data[[j]][[i+2]]
    }
    
  }
  names(covariate.matrices) <- names(covariate.data[[1]])[3:length(covariate.data[[1]])]


#### make matrices for Prim and season
  Prim <- matrix(NA, nrow=sum(num.IDs), ncol=max.time)
  season <- matrix(NA, nrow=sum(num.IDs), ncol=max.time)
  for(j in 1:n.sites){
      Prim[IDsites$start[j]:IDsites$end[j], 1:times[j]] <- matrix(covariate.data[[j]]$temporal.covariates$Prim,nrow=num.IDs[j],ncol=times[j],byrow=TRUE)
      season[IDsites$start[j]:IDsites$end[j], 1:times[j]] <- matrix(covariate.data[[j]]$temporal.covariates$season,nrow=num.IDs[j],ncol=times[j],byrow=TRUE)
    
  }
    
  combined.data <- list.prepend(covariate.matrices, individual.covariates=combined.individual.cov.data )
  combined.data <- list.append(combined.data,Prim=Prim,season=season)
  return(combined.data)
}





combine.monthlyCH.fun <- function(monthlyCH = list(monthlyCH_GC,monthlyCH_N,monthlyCH_Z)){
  
  n.sites <- length(monthlyCH)
  num.IDs <- unlist(lapply(monthlyCH,function(x){dim(x)[1]})) # number of individuals in data sets
  times <- unlist(lapply(monthlyCH,function(x){dim(x)[2]}))
  max.time <- max(times)
  IDsites <- data.frame(site=1:n.sites,start=c(1,rep(NA,n.sites-1)), end=c(num.IDs[1],rep(NA,n.sites-1)))
  for(i in 2:n.sites){
    IDsites$start[i] <- IDsites$end[i-1]+1
    IDsites$end[i] <- IDsites$start[i] + num.IDs[i]-1
  }
  
  new.monthlyCH <- matrix(NA, nrow=sum(num.IDs), ncol=max.time)
    
  for(j in 1:n.sites){
      
    new.monthlyCH[IDsites$start[j]:IDsites$end[j], 1:times[j]] <- monthlyCH[[j]]
    
    
  }
  
  return(new.monthlyCH)
}



combine.porc.fun <- function(p.or.c = list(p.or.c_GC, p.or.c_N, p.or.c_Z)){
  n.sites <- length(p.or.c)
  num.IDs <- unlist(lapply(p.or.c,function(x){dim(x)[1]}))
  months <- unlist(lapply(p.or.c,function(x){dim(x)[2]}))
  days <- unlist(lapply(p.or.c,function(x){dim(x)[3]}))
  IDsites <- data.frame(site=1:n.sites,start=c(1,rep(NA,n.sites-1)), end=c(num.IDs[1],rep(NA,n.sites-1)))
  for(i in 2:n.sites){
    IDsites$start[i] <- IDsites$end[i-1]+1
    IDsites$end[i] <- IDsites$start[i] + num.IDs[i]-1
  }
  
  new.p.or.c <- array(NA, dim=c(sum(num.IDs), max(months), max(days)))
  
  for(j in 1:n.sites){
    
    new.p.or.c[IDsites$start[j]:IDsites$end[j], 1:months[j], 1:days[j]] <- p.or.c[[j]]
    
  }
  
  return(new.p.or.c)
  
}


combine.n.sec.occ.func <- function(n.sec.occ = list(n.sec.occ_GC,n.sec.occ_N,n.sec.occ_Z)){
  n.sites <- length(n.sec.occ)
  num.IDs <- unlist(lapply(n.sec.occ,function(x){dim(x)[1]}))
  months <- unlist(lapply(n.sec.occ,function(x){dim(x)[2]}))
  IDsites <- data.frame(site=1:n.sites,start=c(1,rep(NA,n.sites-1)), end=c(num.IDs[1],rep(NA,n.sites-1)))
  for(i in 2:n.sites){
    IDsites$start[i] <- IDsites$end[i-1]+1
    IDsites$end[i] <- IDsites$start[i] + num.IDs[i]-1
  }
  
  
  new.n.sec.occ <- matrix(NA, nrow=sum(num.IDs), ncol= max(months))
  
  for(j in 1:n.sites){
    
    new.n.sec.occ[IDsites$start[j]:IDsites$end[j], 1:months[j]] <- n.sec.occ[[j]]
    
  }
  
  return(new.n.sec.occ)
  
}



combine.months.trapped.func <- function(months.trapped.mat = list(months.trapped.mat_GC,months.trapped.mat_N,months.trapped.mat_Z)){
  n.sites <- length(months.trapped.mat)
  num.IDs <- unlist(lapply(months.trapped.mat,function(x){dim(x)[1]}))
  months <- unlist(lapply(months.trapped.mat,function(x){dim(x)[2]}))
  IDsites <- data.frame(site=1:n.sites,start=c(1,rep(NA,n.sites-1)), end=c(num.IDs[1],rep(NA,n.sites-1)))
  for(i in 2:n.sites){
    IDsites$start[i] <- IDsites$end[i-1]+1
    IDsites$end[i] <- IDsites$start[i] + num.IDs[i]-1
  }
  
  
  new.months.trapped.mat <- matrix(NA, nrow=sum(num.IDs), ncol= max(months))
  
  for(j in 1:n.sites){
    
    new.months.trapped.mat[IDsites$start[j]:IDsites$end[j], 1:months[j]] <- months.trapped.mat[[j]]
    
  }
  
  return(new.months.trapped.mat)
  
}


cjs.init.z.combined <- function(ch, f) {
  for (i in 1:dim(ch)[1]) {
    max.month <- max(which(!is.na(ch[i,])))
    if (sum(ch[i, 1:max.month]) == 1) next
    n2 <- max(which(ch[i,1:max.month ] == 1))
    ch[i, f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]) {
    ch[i, 1:f[i]] <- NA
  }
  return(ch)
}






####################################################################################################
## Unused 
####################################################################################################

### don't use this one (it's in the combine.covariates.fun)
# combine individual covariates data frames
ind.covariates.paste.fun <- function(individual.covariates = list(covariate.data_GC$individual.covariates,covariate.data_N$individual.covariates, covariate.data_Z$individual.covariates)){
  
  n.sites <- length(individual.covariates)
  site.names <- unlist(lapply(individual.covariates,function(x){unlist(strsplit(as.character(x$tag[1]),"[.]"))[1]}))
  num.IDs <- unlist(lapply(individual.covariates,function(x){dim(x)[1]})) # number of individuals in data sets
  
  # add site covariate
  for(i in 1:n.sites){
    individual.covariates[[i]]$site <- site.names[i]
  }
  
  #change ID so sequential
  for(i in 2:n.sites){
    s <- individual.covariates[[i-1]]$ID[length(individual.covariates[[i-1]]$ID)]+1
    individual.covariates[[i]]$ID <- s:(s+num.IDs[i]-1)
  }
  
  combined.data <- individual.covariates[[1]]
  for(i in 2:n.sites){
    combined.data <- rbind(combined.data,individual.covariates[[i]])
  }
  
  return(combined.data)
}
