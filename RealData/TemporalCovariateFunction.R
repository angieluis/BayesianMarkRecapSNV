#### add months not trapped so can have temporal covariates for all months, and line them up to sessions. so there will be multiple temporal covariates assigned to one primary session if months not trapped.

# This applies to survival (phi) 
# if want to apply covariates to capture probabilities, then months not trapped won't matter, only use the first one


temporaldata.fun <-function(data,site,web){
  sessions <- sort(unique(data$Session[which(data$site==site&data$web==web)]))
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
    
  temporal.covariates <- data.frame(month.session=1:length(ms),session= s1,year=ys,month=ms,covariate.prim=session.num)

  return(temporal.covariates)
}

#temporal.covariates <- temporaldata.fun(data=UNMcaptures, site="Zuni",web=2)




