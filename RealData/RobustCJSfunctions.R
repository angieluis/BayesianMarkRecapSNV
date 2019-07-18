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


## function to make a primary Ch but by week instead of month
# make weeks not trapped = 0 (so not really capture history)
# this will just be used to set initial values and known states for weekly z
# not used:
weekly.primaryCH.fun <- function(CH.primary,temporal.covariates){
  CH <- matrix(NA, ncol=max(temporal.covariates$week),nrow=dim(CH.primary)[1])
  old <- which(is.finite(temporal.covariates$Prim))
  CH[,old]<- CH.primary
  CH <- replace(CH,is.na(CH),0)
  return(CH)
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


### not using:
monthly.temporaldata.fun.old <-function(data,site,web=NULL){
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


library(tidyverse,lubridate) 

# not using:
monthly.temporaldata.fun <-function(
  #dates, # dates trapped
  capture.data, # capture data frame like UNMdata
  temporal.data=NULL, # data frame of monthlytemporal data, with either a column called date or yearmon (must include months not trapped)
  longdata=TRUE, # if TRUE, then data are long like sw.temp.data and have all the sites and all temporal data with no time lags, if false, then data are just as will be input with 
  site=NULL, # if longdata=TRUE 
  web=NULL, # if longdata=TRUE
  cov.list=NULL, # if longdata=TRUE, list of covariates and their time lags, e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag and with a lag 1 and tmax with lag 3. 
  individual.covariates=NULL # if present, will add lists of temporal covariates organized by individual (web) so can be input directly in model, e.g. NDVI_1[i,t]
){
  
  s <- character() # assuming length of site==1 here but not below
  for(i in 1:length(web)){
    s <- c(s,sort(unique(capture.data$Session[which(capture.data$site==site & capture.data$web==web[i])])))
  }
  sessions.trapped <- sort(unique(s))
  
  
  first.session <- sessions.trapped[1]
  last.session <- sessions.trapped[length(sessions.trapped)]
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
  #yearmonth <- paste(as.character(ms),as.character(ys),sep="")
  mc <- as.character(ms)
  for(i in 1:length(ms)){
    mc[i] <- ifelse(ms[i]<10,paste(paste(as.character(ys[i]),"0",as.character(ms[i]),sep="")),paste(paste(as.character(ys[i]),as.character(ms[i]),sep="")))
  }
  
  s1 <- sessions.trapped[match(mc,sessions.trapped)]
  not.trapped <- which(is.na(s1))
  sn <- 1:length(sessions.trapped)
  session.num <- sn[match(mc,sessions.trapped)]
  for(i in 1:length(not.trapped)){
    session.num[not.trapped[i]]<-session.num[max(which(sn<not.trapped[i]))]
  }
  Prim=session.num
  Prim[which(is.na(s1))]<- NA
  
  month.data <- data.frame(long.month=1:length(ms),session= s1,year=ys,month=ms,covariate.prim=session.num,Prim=Prim)
  
  if(length(temporal.data)>0){
    if(longdata==TRUE){
      ls <- length(site)
      datas <- temporal.data[grep(site[1],temporal.data$site,ignore.case=TRUE),]
      if(ls>1){
        for(i in 2:ls){
          datas <- rbind(datas,temporal.data[grep(site[i],temporal.data$site,ignore.case=TRUE),])
        }
      }
      lw <- length(web)
      dataw <- datas[which(datas$web==web[1]),]
      if(lw>1){
        for(i in 2:lw){
          dataw <- rbind(dataw,datas[which(datas$web==web[i]),])
        }
      }
      
      # make a wide data frame with date/yearmon going from first trapped session to last trapped session
      wdate=lubridate::dmy(paste("1",month.data$month,month.data$year,sep="-"))
      
      data.w <- data.frame(date=sort(unique(wdate)))
      data.w$year <- lubridate::year(data.w$date)
      data.w$month <- lubridate::month(data.w$date)
      dataw$date <- lubridate::dmy(paste("1",dataw$yearmon))
      cl <- length(cov.list)
      for(c in 1:cl){
        nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
        col <- grep(names(cov.list)[c],names(dataw),ignore.case=TRUE)
        fd <- which(dataw$date==data.w$date[1]) 
        ld <- which(dataw$date==data.w$date[length(data.w$date)])
        if(length(fd)==1){
          data.w <- cbind(data.w, dataw[(fd-cov.list[[c]]):(ld-cov.list[[c]]),col])
          names(data.w)[dim(data.w)[2]] <- nam
        }
        if(length(fd)>1){
          for(i in 1:length(fd)){
            data.w <- cbind(data.w, dataw[(fd[i]-cov.list[[c]]):(ld[i]-cov.list[[c]]),col])
            names(data.w)[dim(data.w)[2]] <- paste(nam,paste("web",dataw$web[fd[i]],sep=""),sep=".") 
          } 
        }
      }
      month.data <- dplyr::left_join(month.data,data.w) 
    }
  }
  monthly.temporal.data <- list(monthly.longdata=month.data)
  
  # make matrices of temporal covariates matched up to individuals (temporal data for the individual based on which web they were on)
  if(length(individual.covariates)>0){
    for(c in 1:length(cov.list)){
      nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
      cols <- grep(nam,names(month.data))
      col.name <- names(month.data)[cols]
      web.nam <- unlist(lapply(strsplit(col.name,".web"),function(x){x[2]}))
      dat <- month.data[,cols]
      names(dat) <- web.nam 
      
      mat <- matrix(NA,ncol=dim(month.data)[1], nrow=dim(individual.covariates)[1])
      for(i in 1:dim(individual.covariates)[1]){
        w <- individual.covariates$web[i]
        mat[i,] <- dat[,which(names(dat)==w)]
      } #i
      
      
      monthly.temporal.data[[c+1]] <- mat
      names(monthly.temporal.data)[[c+1]] <- nam
    } #c
  } #if
  return(monthly.temporal.data)
}


# creates monthly capture history (including months not trapped) to pass to the
# initial values and known state functions
# make months not trapped = 0 (so not really capture history)
monthly.primaryCH.fun <- function(CH.primary,temporal.covariates){
  CH <- matrix(NA, ncol=max(temporal.covariates$long.month),nrow=dim(CH.primary)[1])
  old <- which(is.finite(temporal.covariates$Prim))
  CH[,old]<- CH.primary
  CH <- replace(CH,is.na(CH),0)
  return(CH)
  
}



####### temporal data on a weekly scale



weekly.temporaldata.fun <-function(
    dates, # dates trapped
    data=NULL, # data frame of monthlytemporal data, with either a column called date or yearmon (must include months not trapped)
    longdata=TRUE, # if TRUE, then data are long like sw.temp.data and have all the sites and all temporal data with no time lags, if false, then data are just as will be input with 
    site=NULL, # if longdata=TRUE 
    web=NULL, # if longdata=TRUE
    cov.list=NULL, # if longdata=TRUE, list of covariates and their time lags, e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag and with a lag 1 and tmax with lag 3. 
    individual.covariates=NULL # if present, will add lists of temporal covariates organized by individual (web) so can be input directly in model, e.g. NDVI_1[i,t]
    ){ 

  time.int <- diff(dates)
  first.dates <- dates[c(1,1+which(time.int>1))]
  
  session.week.data <- data.frame(Prim=1:length(first.dates),week=lubridate::week(first.dates),month=lubridate::month(first.dates),year=lubridate::year(first.dates))
  session.week.data$cumweek<-session.week.data$week+(session.week.data$year-min(session.week.data$year))*52

  time.intervals <- diff(session.week.data$cumweek)

  longweek <- session.week.data$cumweek[1]:session.week.data$cumweek[length(session.week.data$cumweek)]
  longdates <- first.dates[1]
  for(i in 1:length(time.intervals)){
  
    for(j in 1:time.intervals[i]){
      longdates <- c(longdates,longdates[length(longdates)]+lubridate::dweeks(1))
    }
  }

  week.data <- data.frame(week=1:length(longdates),cumweek=longweek,longdates=longdates, month=lubridate::month(longdates),year=lubridate::year(longdates),session=ifelse(lubridate::month(longdates)>9,paste(lubridate::year(longdates),lubridate::month(longdates),sep=""),paste(lubridate::year(longdates),"0",lubridate::month(longdates),sep="")))
  week.data$Prim=session.week.data$Prim[match(week.data$cumweek,session.week.data$cumweek)]
  
  if(length(data)>0){
  if(longdata==TRUE){
    ls <- length(site)
    datas <- data[grep(site[1],data$site,ignore.case=TRUE),]
    if(ls>1){
      for(i in 2:ls){
        datas <- rbind(datas,data[grep(site[i],data$site,ignore.case=TRUE),])
       }
     }
    lw <- length(web)
    dataw <- datas[which(datas$web==web[1]),]
    if(lw>1){
      for(i in 2:lw){
        dataw <- rbind(dataw,datas[which(datas$web==web[i]),])
      }
    }
    
    # make a wide data frame with date/yearmon going from first trapped session to last trapped session
    wdate=lubridate::dmy(paste("1",week.data$month,week.data$year,sep="-"))
    data.w <- data.frame(date=sort(unique(wdate)))
    data.w$year <- lubridate::year(data.w$date)
    data.w$month <- lubridate::month(data.w$date)
    dataw$date <- lubridate::dmy(paste("1",dataw$yearmon))
    cl <- length(cov.list)
    for(c in 1:cl){
      nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
      col <- grep(names(cov.list)[c],names(dataw),ignore.case=TRUE)
      fd <- which(dataw$date==data.w$date[1]) 
      ld <- which(dataw$date==data.w$date[length(data.w$date)])
      if(length(fd)==1){
        data.w <- cbind(data.w, dataw[(fd-cov.list[[c]]):(ld-cov.list[[c]]),col])
        names(data.w)[dim(data.w)[2]] <- nam
      }
      if(length(fd)>1){
        for(i in 1:length(fd)){
          data.w <- cbind(data.w, dataw[(fd[i]-cov.list[[c]]):(ld[i]-cov.list[[c]]),col])
          names(data.w)[dim(data.w)[2]] <- paste(nam,paste("web",dataw$web[fd[i]],sep=""),sep=".") 
        } 
      }
    }
  
    week.data <- dplyr::left_join(week.data,data.w)  
    
  }

  }
  
  weekly.temporal.data <- list(weekly.longdata=week.data)

  # make matrices of temporal covariates matched up to individuals (temporal data for the individual based on which web they were on)
  if(length(individual.covariates)>0){
    for(c in 1:length(cov.list)){
      nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
      cols <- grep(nam,names(week.data))
      col.name <- names(week.data)[cols]
      web.nam <- unlist(lapply(strsplit(col.name,".web"),function(x){x[2]}))
      dat <- week.data[,cols]
      names(dat) <- web.nam 
      
      mat <- matrix(NA,ncol=dim(week.data)[1], nrow=dim(individual.covariates)[1])
      for(i in 1:dim(individual.covariates)[1]){
        w <- individual.covariates$web[i]
        mat[i,] <- dat[,which(names(dat)==w)]
      } #i
      
       
      weekly.temporal.data[[c+1]] <- mat
      names(weekly.temporal.data)[[c+1]] <- nam
    } #c
  } #if
  return(weekly.temporal.data)
}





###############################
## Create a dataframe for individual covariates
# this has 2 individual covariates: web and sex
# tags must line up to Capture history rows
individual.covariate.fun <- function(data, #like UNMdata
                                     tags  # must line up w/ CH.secondary rows
                                     ){ 
  nind <- length(tags)
  
  ic <- data.frame(ID=1:nind,tag=tags)
  
  web <- character()
  sex <- numeric() #1 male, 0 female
  for(i in 1:nind){
    ind <- which(data$tag==IDs[i])
    x <- data[ind,]
    x <- x[order(x$Session),]
    web[i] <- as.character(x$web[1])
    
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- factor(web)
  ic$sex <- factor(replace(sex,sex==-9,NA))
  return(ic)
}



######################################################
# function to take Robust Design capture history list of primary 
# occasions (months) and turn it into long data frame where each 
# row is a month (including months not trapped)

monthly.longdataCH.fun<-function(CH.secondary, 
                                 temporal.covariate.df=NULL, 
                                 individual.covariates=NULL, 
                                 p_or_c=FALSE){ 
  #  add the primary occassion to the data (if not already there)
  if(length(which(colnames(CH.secondary[[1]])=="Prim"))==0){
    for(i in 1:length(CH.secondary)){
      CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
    }
  }
  #### define first capture (f) if not already there
  # z will include months not trapped, so first.caught needs to reflect long.month not Prim
  if(length(which(names(individual.covariates)=="f.longmonth"))==0){
    CH.primary <- primary.ch.fun(CH.secondary)
    first.caught <- apply(CH.primary,1,function(x){min(which(x>0))}) #gives primary occasion first caught (not long.month)
    individual.covariates$f.longmonth <- temporal.covariate.df$long.month[match(first.caught,temporal.covariate.df$Prim)]
  }

  obs.dat <- purrr::map_df(
    CH.secondary, 
    ~ tibble::as_tibble(.x) %>%    # 
      dplyr::mutate(
        ID = 1:n()
      ) %>% 
      tidyr::gather(Sec, State, -Prim, -ID) %>%
      dplyr::select(ID, Prim, Sec, State)
  )
  

  # join with temporal covariates to make monthly
  obs.dat.full <- inner_join(obs.dat,temporal.covariate.df[,c("long.month","Prim")])
  obs.dat.full <- arrange(obs.dat.full,long.month,ID)
  
  #  Subset observation data to observed bits
  obs.dat.full <- dplyr::left_join(obs.dat.full, individual.covariates[,c("ID","f.longmonth")]) %>%
    dplyr::filter(long.month >= f.longmonth)
  
  return(obs.dat.full)
  
}


######################################################
# function to take Robust Design capture history list of primary 
# occasions (months) and turn it into long data frame where each 
# row is a week

weekly.longdataCH.fun<-function(CH.secondary, temporal.covariates, individual.covariates, p_or_c=FALSE){ 
  #  add the primary occassion to the data
  for(i in 1:length(CH.secondary)){
    CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
  }

  #### define first capture (f)
  # z is going to be by week - but capture histories are by month. need f to be in weeks
  first.caught <- apply(CH.primary,1,function(x){min(which(x>0))})
  individual.covariates$f.week <- temporal.covariates$week[match(first.caught,temporal.covariates$Prim)]

  obs.dat <- purrr::map_df(
    CH.secondary, 
    ~ tibble::as_tibble(.x) %>%    # 
      dplyr::mutate(
        ID = 1:n()
      ) %>% 
      tidyr::gather(Sec, State, -Prim, -ID) %>%
      dplyr::select(ID, Prim, Sec, State)
  )

  #### if p_or_c==TRUE, add a column that indicates whether that 
  #individual has been caught before in that primary occasion (do we 
  # use p or c?)
  if(p_or_c==TRUE){
    p.or.c <- numeric()
    for(i in 1:dim(obs.dat)[1]){
      # the times that animal was caught that primary session
      dat <- obs.dat[which(obs.dat$Prim==obs.dat$Prim[i] & obs.dat$ID==obs.dat$ID[i] & obs.dat$State==1),]
  
      if(dim(dat)[1]==0){ # if not caught that primary session at all use p (0)
        p.or.c[i] <- 0
      }else{ #  otherwise use p (0) unless already caught caught that session, then use c (1)
        firstcap <- min(as.numeric(dat$Sec))
        p.or.c[i] <- ifelse(firstcap<obs.dat$Sec[i], 1 ,0)   #BUGs doesn't like characters so 0 is p, 1 is c
      }
    }
    p.or.c <- factor(p.or.c)
    obs.dat <- dplyr::mutate(obs.dat, p.or.c)
  }

  # join with temporal covariates to make weekly
  obs.dat.full <- inner_join(obs.dat,temporal.covariates[,c("week","Prim")])
  obs.dat.full <- arrange(obs.dat.full,week,ID)

  #  Subset observation data to observed bits
  obs.dat.full <- dplyr::left_join(obs.dat.full, individual.covariates[,c("ID","f.week")]) %>%
    dplyr::filter(week >= f.week)

  obs.dat.full$Sec=as.numeric(obs.dat.full$Sec)
  return(obs.dat.full)
}



######################################
## functions for plotting output

hist.plot.fun <- function(BUGSout, ...){
  par(mfrow=c(1,1))
  hist(BUGSout$BUGSoutput$sims.list[[1]],breaks = 30,main=names(BUGSout$BUGSoutput$sims.list)[[1]],xlab="value")
  par(ask=TRUE)
  for(i in 2:length(names(BUGSout$BUGSoutput$sims.list))){
    if(dim(BUGSout$BUGSoutput$sims.list[[i]])[2]==1){
      hist(BUGSout$BUGSoutput$sims.list[[i]],breaks = 30,main=names(BUGSout$BUGSoutput$sims.list)[[i]],xlab="value")
    }
    if(dim(BUGSout$BUGSoutput$sims.list[[i]])[2]>1){
      for(j in 1:dim(BUGSout$BUGSoutput$sims.list[[i]])[2] ){
        hist(BUGSout$BUGSoutput$sims.list[[i]][,j],breaks = 30,main=paste(names(BUGSout$BUGSoutput$sims.list)[[i]],j),xlab="value")
      }
    }
  }
  par(ask=FALSE)
}

chain.plot.fun <- function(BUGSout, ...){
  pari <- (1:length(BUGSout$BUGSoutput$sims.list))[-which(names(BUGSout$BUGSoutput$sims.list)=="deviance")]
  par(mfrow=c(1,1))
  plot.ts(BUGSout$BUGSoutput$sims.list[[1]],main=names(BUGSout$BUGSoutput$sims.list)[[1]],ylab="value")
  par(ask=TRUE)
  for(i in 2:length(pari)){
    if(dim(BUGSout$BUGSoutput$sims.list[[i]])[2]==1){
      plot.ts(BUGSout$BUGSoutput$sims.list[[i]],main=names(BUGSout$BUGSoutput$sims.list)[[i]],ylab="value")
    }
    if(dim(BUGSout$BUGSoutput$sims.list[[i]])[2]>1){
      for(j in 1:dim(BUGSout$BUGSoutput$sims.list[[i]])[2] ){
        plot.ts(BUGSout$BUGSoutput$sims.list[[i]][,j],main=paste(names(BUGSout$BUGSoutput$sims.list)[[i]],j),ylab="value")
      }
    }
  }
  par(ask=FALSE)
}





# make p.or.c array, with a 0 if not caught yet that session and 1 if have
p.or.c.array.fun<- function(CH.secondary, # can be list or array
                            temporal.covariates, #dataframe from list
                            n.sec.occ, # number of secondary occasions
                            list=FALSE, #is it a list or array?
                            weekly=FALSE
                            ){ #is it weekly or monthly?
  nind <- ifelse(list==TRUE,dim(CH.secondary[[1]])[1],dim(CH.secondary)[1])
  
  temp.data <- temporal.covariates
  
  nt <- dim(temp.data)[1]
  n.sec <- n.sec.occ    
  p.or.c <- array(NA, dim=c(nind, nt, max(n.sec)))
  
  if(weekly==TRUE){
   weeks.trapped <- temp.data$week[which(is.finite(temp.data$Prim))]
    for(w in weeks.trapped){
      m <- temp.data$Prim[w]
      for(i in 1:nind){
        for(d in 1:n.sec[m]){
          p.or.c[i,w,d] <- ifelse(sum(CH.secondary[[m]][i,1:(d-1)])==0,0,1)
        } #d
      } #i
    } #w
  } else{
    t.trapped <- temp.data$long.month[which(is.finite(temp.data$Prim))]
    for(t in t.trapped){
      m <- temp.data$Prim[t]
      for(i in 1:nind){
        for(d in 1:n.sec[m]){
          dsum <- ifelse(list==TRUE,sum(CH.secondary[[m]][i,1:(d-1)]), sum(CH.secondary[i,m,1:(d-1)]))
          p.or.c[i,t,d] <- ifelse(dsum==0,0,1)
        } #d
      } #i
    } #w
  }
  
  return(p.or.c)
}


#############################################################
## all in one function to return all covariate data
## output is a list with elements:
# [[1]]$individual.covariates  #with tag, web, sex, month first caught
# [[2]]$temporal.covariates # with session, long.month, month covariates, etc
# after that are matrices of dimensions [i,m] that have the temporal covariates
# aligned to individuals and long.month (all months not just those trapped)
# these matrices can be plugged into models for phi[i,m]



library(tidyverse,lubridate,stringr) 
#source('~/Documents/JAGS/BayesianMarkRecapSNV/RealData/DiversityFunctions.R')
# might want to put those functions in here too

monthly.covariate.function <-function(
  capture.data, # capture data frame like Zuni12.pema.data. best if cut down 
                # to data used so no overlap with tags among sites/species
                # or specify tags by site.tag and say so below
  CH.secondary, # as monthly list
  tags, # tag names that line up to CH.secondary
  by.sitetag=FALSE, # if TRUE, then tags are specified by site e.g., "Zuni.1101"
  sessions, # all sessions to include e.g. 199806 (sessions trapped even if no pm caught)
  temporal.data=NULL, # data frame of monthly temporal data, with either a column 
                      # called date or yearmon (must include months not trapped), 
                      # must be in long format like sw.temp.data and have all the
                      # sites and all temporal data with no time lags
  diversity.data=NULL, # data frame of species MNAs and diversities from 
                      # diversity.df.function(), species 2 letter codes in lowercase
                      # cov.list below must match column headers
                      # input should be scaled (scale=TRUE in fxn)
  site=NULL, #  e.g. "Zuni"
  web=NULL, #  e.g. c("1","2"), 
  cov.list=NULL # list of temporal covariates and their time lags, 
                # e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag 
                # and with a lag 1 and tmax with lag 3. 
){
  
  names(capture.data) <- tolower(names(capture.data))
  nind <- length(tags)
  #site.tag <- paste(str_to_lower(site),tags,sep=".")
  ic <- data.frame(ID=1:nind,tag=tags)
  
  webi <- character()
  sex <- numeric() #1 male, 0 female
  for(i in 1:nind){
    if(by.sitetag==FALSE){
      ind <- which(capture.data$tag==tags[i])
    } else{
      tagi <- str_split(tags[i],"[.]")[[1]][2]
      sitei <- str_split(tags[i],"[.]")[[1]][1]
      ind <- which(capture.data$tag==tagi & str_to_lower(capture.data$site)==sitei)
    }
    x <- capture.data[ind,]
    x <- x[order(x$session),]
    webi[i] <- as.character(x$web[1])
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- factor(webi)
  ic$sex <- replace(sex,sex==-9,0.5) #split the difference for unknown sexes
  
  sessions.trapped <- sort(unique(sessions))
  
  first.session <- sessions.trapped[1]
  last.session <- sessions.trapped[length(sessions.trapped)]
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
  #yearmonth <- paste(as.character(ms),as.character(ys),sep="")
  mc <- as.character(ms)
  for(i in 1:length(ms)){
    mc[i] <- ifelse(ms[i]<10,paste(paste(as.character(ys[i]),"0",as.character(ms[i]),sep="")),paste(paste(as.character(ys[i]),as.character(ms[i]),sep="")))
  }
  
  s1 <- sessions.trapped[match(mc,sessions.trapped)]
  not.trapped <- which(is.na(s1))
  sn <- 1:length(sessions.trapped)
  session.num <- sn[match(mc,sessions.trapped)]
  for(i in 1:length(not.trapped)){
    session.num[not.trapped[i]]<-session.num[max(which(sn<not.trapped[i]))]
  }
  Prim=session.num
  Prim[which(is.na(s1))]<- NA
  
  month.data <- data.frame(long.month=1:length(ms),session= s1,year=ys,month=ms,covariate.prim=session.num,Prim=Prim)
  
  if(length(temporal.data)>0){ # |length(diversity.data)>0
    
      ls <- length(site)
      datas <- temporal.data[grep(site[1],temporal.data$site,ignore.case=TRUE),]
      if(ls>1){
        for(i in 2:ls){
          datas <- rbind(datas,temporal.data[grep(site[i],temporal.data$site,ignore.case=TRUE),])
        }
      }
      lw <- length(web)
      dataw <- datas[which(datas$web==web[1]),]
      if(lw>1){
        for(i in 2:lw){
          dataw <- rbind(dataw,datas[which(datas$web==web[i]),])
        }
      }
      
      
      # make a wide data frame with date/yearmon going from first trapped session to last trapped session
      wdate=lubridate::dmy(paste("1",month.data$month,month.data$year,sep="-"))
      
      data.w <- data.frame(date=sort(unique(wdate)))
      data.w$year <- lubridate::year(data.w$date)
      data.w$month <- lubridate::month(data.w$date)
      dataw$date <- lubridate::dmy(paste("1",dataw$yearmon))
      
      # now paste in diversity data if present 
      if(length(diversity.data)>0){
        #use dates to line up diversity data to temporal data, then add NAs elswhere
        diversity.data$date <- lubridate::dmy(paste("1",diversity.data$month,diversity.data$year,sep="-"))
        diversity.data$site <- str_to_lower(diversity.data$site)
        dataw <- dplyr::left_join(dataw,diversity.data)
      }
      
      
      cl <- length(cov.list)
      for(c in 1:cl){
        nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
        col <- grep(names(cov.list)[c],names(dataw),ignore.case=TRUE)
        fd <- which(dataw$date==data.w$date[1]) 
        ld <- which(dataw$date==data.w$date[length(data.w$date)])
        if(length(fd)==1){
          data.w <- cbind(data.w, dataw[(fd-cov.list[[c]]):(ld-cov.list[[c]]),col])
          names(data.w)[dim(data.w)[2]] <- nam
        }
        if(length(fd)>1){
          for(i in 1:length(fd)){
            data.w <- cbind(data.w, dataw[(fd[i]-cov.list[[c]]):(ld[i]-cov.list[[c]]),col])
            names(data.w)[dim(data.w)[2]] <- paste(nam,paste("web",dataw$web[fd[i]],sep=""),sep=".") 
          } 
        }
      }
      month.data <- dplyr::left_join(month.data,data.w)
      
      ind.na <- which(is.na(month.data),arr.ind = TRUE)
      ind.na <- ind.na[-which(ind.na[,2]<8),] # these columns have NAs that need to be filled for CJS models to run. They are likely the first time points when there is a lag. So fill in the next value.
      
      if(dim(ind.na)[1]>0){
        for(i in 1:dim(ind.na)[1]){
          #find the next finite value after this one and plug it in
          vi <- which(is.finite(month.data[,ind.na[i,2]]))
          vin <- vi[min(which(vi > ind.na[i,1]))]
          month.data[ind.na[i,1],ind.na[i,2]] <- month.data[vin,ind.na[i,2]]
        }
      }
    }
  #}
  
  ## add first capture to individual covariates data
  CH.primary <- primary.ch.fun(CH.secondary)
  first.caught <- apply(CH.primary,1,function(x){min(which(x>0))}) #gives primary occasion first caught (not long.month)
  ic$f.longmonth <- month.data$long.month[match(first.caught,month.data$Prim)]
  
  individual.covariates=ic
  covariate.data <- list(individual.covariates=individual.covariates, temporal.covariates=month.data)
  
  # make matrices of temporal covariates matched up to individuals (temporal data for the individual based on which web they were on): dimensions [i,m], where m is longmonth - all months not just those trapped
  
  for(c in 1:length(cov.list)){
    nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
    cols <- grep(nam,names(month.data))
    col.name <- names(month.data)[cols]
    web.nam <- unlist(lapply(strsplit(col.name,".web"),function(x){x[2]}))
    dat <- month.data[,cols]
    names(dat) <- web.nam 
    
    mat <- matrix(NA,ncol=dim(month.data)[1], nrow=dim(individual.covariates)[1])
    for(i in 1:dim(individual.covariates)[1]){
      w <- individual.covariates$web[i]
      mat[i,] <- dat[,which(names(dat)==w)]
    } #i
    
    
    covariate.data[[c+2]] <- mat
    names(covariate.data)[[c+2]] <- nam
  } #c
  
  return(covariate.data)
}



######### function to put NAs into CH.secondary for individuals on webs not trapped that session
# e.g. month 4 wasn't trapped at Zuni2 but was at Zuni2 (or at least no animals were caught)
# use individual covariates to see which individuals/rows were on webs not trapped

# not done, need to test
CH.addNA.function <- function(
         CH.secondary, # list of sessions
         individual.covariates, # dataframe from monthly.covariate.function   $individual.covariates
         webmonths # list with each element a vector of months trapped, starting from 1
         ){
 
  nwebs <- length(webmonths)
  allmonths <- sort(unique(unlist(webmonths)))
  
  for(i in 1:nwebs){
    web <- strsplit(names(webmonths)[i],".")[2]
    not.trapped <- allmonths[which(is.na(match(allmonths,webmonths[[i]])))]
    id.web <- individual.covariates$tag[which(individual.covariates$web==web)]
    inds <- match(id.web,rownames(CH.secondary[[1]]))
    for(m in 1:length(not.trapped)){
      CH.secondary[[not.trapped[m]]][inds,] <- rep(NA,dim(CH.secondary[[not.trapped[m]]])[2])
    }
  }
  
  return(CH.secondary)
  
  
}