
# make sure the temp data is normalized.

multisite.monthly.covariate.function <-function(
  cleaned.data = pema.final.clean, 
  CH.secondary = Ch.list, # as monthly list
  tags = rownames(Ch.list[[1]]), # tag names that line up to CH.secondary
  by.sitetag=TRUE, # if TRUE, then tags are specified by site e.g., "Zuni.1101"
  sessions = session.list$all.sessions, # all sessions to include e.g. 199806 (sessions trapped even if no pm caught)
  temporal.data= sw.temp.data, # data frame of monthly temporal data, with either a column 
  # called date or yearmon (must include months not trapped), 
  # must be in long format like sw.temp.data and have all the
  # sites and all temporal data with no time lags
  diversity.data=NULL, # data frame of species MNAs and diversities from 
  # diversity.df.function(), species 2 letter codes in lowercase
  # cov.list below must match column headers
  # input should be scaled (scale=TRUE in fxn)
  site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"), #  e.g. c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2")
  cov.list=list(ndvi=8,prcp=0,tmax=0) # list of temporal covariates and their time lags, 
  # e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag 
  # and with a lag 1 and tmax with lag 3. 
){
  
  names(cleaned.data) <- tolower(names(cleaned.data))
  tags <- tolower(tags)
  site.webs <- tolower(site.webs)
  cleaned.data$site.tag <- paste(tolower(cleaned.data$site),tolower(cleaned.data$tag),sep=".")
  cleaned.data$site.web <- paste(tolower(cleaned.data$site),tolower(cleaned.data$web),sep=".")
  temporal.data$site.web <- paste(tolower(temporal.data$site),tolower(temporal.data$web),sep=".")
  nind <- length(tags)
  
  ind.clean <- numeric()
  for(i in 1:length(site.webs)){
    ind.clean <- c(ind.clean,which(cleaned.data$site.web==site.webs[i]))
  }
  cleaned.data <- cleaned.data[ind.clean,]
  
  
  ic <- data.frame(ID=1:nind,tag=tags)
  
  site.webi <- character()
  sex <- numeric() #1 male, 0 female
  for(i in 1:nind){
    if(by.sitetag==FALSE){
      ind <- which(cleaned.data$tag==tags[i])
    } else{
      ind <- which(cleaned.data$site.tag==tags[i])
    }
    x <- cleaned.data[ind,]
    x <- x[order(x$session),]
    site.webi[i] <- as.character(x$site.web[1])
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- factor(site.webi)
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
    
    ls <- length(site.webs)
    datas <- temporal.data[which(temporal.data$site.web==site.webs[1]),]
    if(ls>1){
      for(i in 2:ls){
        datas <- rbind(datas,temporal.data[which(temporal.data$site.web==site.webs[i]),])
      }
    }
    
    
    # make a wide data frame with date/yearmon going from first trapped session to last trapped session
    wdate=lubridate::dmy(paste("1",month.data$month,month.data$year,sep="-"))
    
    data.w <- data.frame(date=sort(unique(wdate)))
    data.w$year <- lubridate::year(data.w$date)
    data.w$month <- lubridate::month(data.w$date)
    datas$date <- lubridate::dmy(paste("1",datas$yearmon))
    
    # now paste in diversity data if present 
    if(length(diversity.data)>0){
      #use dates to line up diversity data to temporal data, then add NAs elswhere
      diversity.data$date <- lubridate::dmy(paste("1",diversity.data$month,diversity.data$year,sep="-"))
      diversity.data$site <- str_to_lower(diversity.data$site)
      diversity.data$site.web <- paste(diversity.data$site,tolower(diversity.data$web))
      datas <- dplyr::left_join(datas,diversity.data)
    }
    
    
    cl <- length(cov.list)
    for(c in 1:cl){
      nam <- paste(names(cov.list)[c],cov.list[[c]],sep="_")
      col <- grep(names(cov.list)[c],names(datas),ignore.case=TRUE)
      fd <- which(datas$date==data.w$date[1]) 
      ld <- which(datas$date==data.w$date[length(data.w$date)])
      if(length(fd)==1){
        data.w <- cbind(data.w, datas[(fd-cov.list[[c]]):(ld-cov.list[[c]]),col])
        names(data.w)[dim(data.w)[2]] <- nam
      }
      if(length(fd)>1){
        for(i in 1:length(fd)){
          data.w <- cbind(data.w, datas[(fd[i]-cov.list[[c]]):(ld[i]-cov.list[[c]]),col])
          names(data.w)[dim(data.w)[2]] <- paste(nam,paste("web",datas$site.web[fd[i]],sep=""),sep=".") 
        } 
      }
    }
    month.data <- dplyr::left_join(month.data,data.w)
    
    ind.na <- which(is.na(month.data),arr.ind = TRUE)
    ind.na <- ind.na[-which(ind.na[,2]<8),] # these columns have NAs that need to be filled for CJS models to run. They are likely the first time points when there is a lag. So fill in the next value. (this is after removing the first 7 column because that is just the session numbers etc. not the covariate data)
    
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





# need to make this handle when some webs are trapped different primary and secondary occasions

# make p.or.c array, with a 0 if not caught yet that session and 1 if have
multisite.p.or.c.array.fun<- function(CH.secondary, # can be list or array
                            temporal.covariates #dataframe from list
                            
){ 
  Ch.primary <- primary.ch.fun(CH.secondary)
  nind <- dim(CH.secondary[[1]])[1]
  
  temp.data <- temporal.covariates
  
  nt <- dim(temp.data)[1]
  n.sec <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))
  p.or.c <- array(NA, dim=c(nind, nt, max(n.sec)))
  
  for(i in 1:nind){
    m.trapped <- which(is.finite(Ch.primary[i,])) #sessions this indiv could have been trapped
    longm.trapped <- temp.data$long.month[match(m.trapped,temp.data$Prim)] #what long months
    for(m in 1:length(m.trapped)){
      nsec <- length(which(is.finite(CH.secondary[[m.trapped[m]]][i,])))# number of secondary occasions there was trapping just for that individual(varies by web) that month 
        for(d in 1:nsec){
          dsum <- sum(CH.secondary[[m.trapped[m]]][i,1:(d-1)])
          p.or.c[i,longm.trapped[m],d] <- ifelse(dsum==0,0,1)
        } #d
      } #m
    } #i
  
  
  return(p.or.c)
}


