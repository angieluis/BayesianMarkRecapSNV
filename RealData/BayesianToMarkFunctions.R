

#####################################################################
#####################################################################
#####################################################################

# create a Robust Design capture history to input in Program Mark
# from CH.secondary list used in Bayesian analysis

RDMARKch.from.CHsecondary <- function(CH.secondary){
  
  n.prim <- length(CH.secondary)
  n.ind <- dim(CH.secondary[[1]])[1]
  CH <- matrix(,nrow=nind)
  for(m in 1:n.prim){
    CH <- cbind(CH,CH.secondary[[m]])
  }
  CH <- CH[,-1]
  
  char.ch <- character(length=n.ind)
  for(i in 1:n.ind){
    char.ch[i] <- paste(CH[i,],collapse = "")
  }
  
  RD.ch <- data.frame(ch=char.ch, freq=rep(1,n.ind))
  return(RD.ch)
}


#####################################################################

time.int.from.sessions <- function(CH.secondary,sessions){
  
  n.sec.occ <- unlist(lapply(CH.secondary,function(x){dim(x)[2]}))
  
  year <- numeric()
  month <- numeric()
  for(t in 1:length(sessions)){
    y <- strsplit(as.character(sessions[t]),split=character(0))[[1]][1:4]
    year[t] <- as.numeric(paste(y[1],y[2],y[3],y[4],sep=""))
      
    m <- strsplit(as.character(sessions[t]),split=character(0))[[1]][5:6]
    month[t] <- as.numeric(paste(m[1],m[2],sep=""))
    
  }
  wdate=lubridate::dmy(paste("1",month,year,sep="-"))
  
  month.intervals <- round(diff(wdate)/30)
  
  time.intervals <- rep(0,(n.sec.occ[1]-1))
  for(t in 2:length(n.sec.occ)){
    time.intervals <- c(time.intervals,month.intervals[t-1],
                        rep(0,(n.sec.occ[t]-1)))
  }
  return(time.intervals)
}




#####################################################################
#####################################################################
#####################################################################

## old code for Program Mark straight from original dataset
RDMScapfun=function(data,site,web,species){
  #### function to create Robust Design Multistate (ab neg and pos) capture histories and time intervals from data frame of captures
  ## output is a list with 
  #[[1]] "ch" capture history, where "1" means ab negative and "2" means ab positive
  #[[3]] "uniquedates", all the dates in the dataset
  #[[4]] "time.int" time intervals (months) #this will need to rounded before input in RMark
  
  
  
  data=data[grep(site,data$site,ignore.case=TRUE),]
  data=data[which(data$web==web),]
  date=as.Date(gsub(" ", "",as.character(data$Date)),format="%m/%d/%Y")
  sessions=sort(unique(data$Session))
  uniquedates=sort(unique(date))
  
  data=data[grep(species,data$letter_2,ignore.case=TRUE),]
  
  dt=numeric()
  for(i in 2:length(uniquedates)){
    dt[i]=difftime(uniquedates[i], uniquedates[i-1])
  }
  dt2=dt-2
  dt2=dt2/30.5 # days a month
  dt2=replace(dt2,which(dt2<0),0)
  
  tags=as.character(sort(unique(data$tag)))
  if(length(which(tags=="-9"))>0){
    tags=tags[-which(tags=="-9")]}
  if(length(which(tags==""))>0){
    tags=tags[-which(tags=="")]}
  
  ch=matrix(NA,ncol=length(uniquedates),nrow=length(tags),dimnames=list(tags, uniquedates))
  
  date2=as.Date(gsub(" ", "",as.character(data$Date)),format="%m/%d/%Y")
  
  for(i in 1:length(tags)){
    for(j in 1:length(uniquedates)){
      ch[i,j] =length(which(data$tag==tags[i]&date2== uniquedates[j]))
    }
  }
  
  
  ch=replace(ch,ch>1,1) # this is the regular robust design capture history
  
  #now need to add strata of snv negative and positive
  pos.tags=data$tag[which(data$snv_pos==1)] #which animals were ever positive
  
  for(i in 1:length(pos.tags)){ #loop through each positive animal
    tag.date=sort(date2[data$tag==pos.tags[i]]) #all dates that animal was caught
    first.date.pos=min(date2[which(data$tag==pos.tags[i]&data$snv_pos==1)]) #first date that animal was snv positive
    later.dates=tag.date[which(tag.date>first.date.pos)]# dates after first detected as pos
    # was there a date before in the same session that was negative?
    x=seq(first.date.pos-2,first.date.pos+2,"day")
    y=tag.date[match(x,tag.date)] #which dates were the animal caught within 2 days of being positive
    dates.pos=sort(unique(c(later.dates,y))) #later dates or dates before within the same session should be counted as positive
    
    if(length(dates.pos)>0){ #make those dates positive
      for(j in 1:length(dates.pos)){
        ch[tags==pos.tags[i],which(uniquedates==dates.pos[j])]=2 
      }
    }
  }
  
  # Mark doesn't like it if multiple strata recorded for same animal in same primary session, which could happen if  blood wasn't taken at the first capture for some reason but was at second and tested positive. How many of those are there?
  
  out=list(ch,uniquedates,dt2)
  names(out)=c("ch","uniquedates","time.int")
  return(out)
}



#date()
#RDMSch.Zuni2=RDMScapfun(data=UNMdata,site="Zuni",web=2,species="PM")
#date()


#####################################################################
#####################################################################
#####################################################################


# make MARK design data with all the temporal covariates I have

design.data.function <- function(MARK.process, covariate.data){
  
}
