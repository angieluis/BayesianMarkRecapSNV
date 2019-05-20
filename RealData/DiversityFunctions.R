################################################################
# function to calculate MNAs (in dataframe with sessions & 
# longmonths, including sessions not trapped as NAs)
################################################################

library(zoo)
MNA.function <- function(data=UNMdata, # capture data 
                         site="Zuni", # 1 site at a time
                         web=1, # 1 web at a time
                         species="pm", # can put in more than 1
                         sessions=session.list$web.1 # sessions trapped
                         ){
  # output is dataframe with columns:
  # long.month: from 1 to number of total months spanned (could be more than months trapped)
  # session
  # year
  # month
  # Prim: primary occasion from 1 to number of months trapped
  # MNA: this will include NAs for months not trapped
  # MNA.interp : interpolated values for NAs (mean of previous & next)
  
  data=data[grep(site,data$site,ignore.case=TRUE),]
  data1=data[which(data$web==web),]
  
  data=data1[grep(species[1],data1$letter_2,ignore.case=TRUE),]
  if(length(species)>1){
    for(i in 2:length(species)){
      data=rbind(data,data1[grep(species[i],data1$letter_2,ignore.case=TRUE),])
    }
  }
  
  tags=as.character(sort(unique(data$tag)))
  if(length(which(tags=="-9"))>0){
    tags=tags[-which(tags=="-9")]}
  
  #######  X:  all the code between X's is to create a full length of months & sessions that span the time frame and include those not trapped
  sessions.trapped <- sessions
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
  
  #######  X
  
  
  # since this will be used in a function for diversity, some of the species might not have been trapped or have tag numbers so return 0s
  if(length(tags)==0){
    MNA <- rep(0, length(sessions.trapped))
  } else{
    
  ch=matrix(NA,ncol=length(sessions),nrow=length(tags),dimnames=list(tags,sessions))
  
  
  for(i in 1:length(tags)){
    for(j in 1:length(sessions)){
      ch[i,j] =length(which(data$tag==tags[i]&data$Session==sessions[j]))
    }
  }
  
  #plot.ts(apply(ch,2,sum))
  
  ch2=replace(ch,ch>1,1)
  #lines(apply(ch2,2,sum))
  
  chMNA=ch2
  
  for(i in 1:length(tags)){
    x=which(ch2[i,]>0)
    if(length(x)>1){
      if(min(x)<(max(x)-1)){
        chMNA[i,min(x):max(x)]=1			
      }
    }
  }
  
  MNA=apply(chMNA,2,sum)

  }
  

  MNA.data <- data.frame(long.month=1:length(ms),session= s1,year=ys,month=ms,Prim=Prim, MNA=MNA[match(s1,sessions)], MNA.interp = na.approx(MNA[match(s1,sessions)],1:length(s1)))
  
  return(MNA.data)
  
}





################################################################
# function to create dataframe of all species MNAs for a site
# and diversity indices
### interpolating the NAs
## need to make this work for multiple webs at a site
# need to make the interpolation optional
################################################################
#Using the vegan package to calculate diversity indices
library(vegan)
library(stringr)
library(dplyr)
diversity.df.function <- function(
        data=UNMdata, # capture data 
        sites="Zuni", # 1 site at a time
        webs=1, # trying for multiple webs
        sessions=session.list$web.1, # sessions trapped - may want to change to specify different sessions for different webs? (so will put in NA before or after ever trapping?)
        interpolate=FALSE, # do you want to interpolate NAs?
        scale=FALSE, # do you want to scale each one between 0 and 1 (divide by max) - this is needed for CJS models
        include.pm = TRUE # reviewers wanted diversity indices calculated with pm removed (include.pm=FALSE), but removing pm leads to Inf in invSimpsonD calculations (when no other species are present)
                                  
){
  
  data <- filter(data, site %in% sites, web %in% webs )

  species <- sort(unique(str_to_lower(data$letter_2)))
  species <- species[-which(species=="-9")]
  species <- species[-which(species=="")]
  # check that each of these species has at least one tag number associated with it
  t <- numeric()
  for(i in 1:length(species)){
    dt <- sort(unique(data$tag[grep(species[i],data$letter_2,ignore.case=TRUE)]))
    baddies <- which(dt=="-9"|dt==""|is.na(dt))
    if(length(baddies)>0){
      dt <- dt[-baddies]
    }
    t[i] <- length(dt)
  }
  
  #remove species with no tag numbers
  sb <- which(t==0) 
  if(length(sb)>0){
    species <- species[-sb]
  }   
  
  
  
  w <- 1
  while(w <= length(webs)){
    # set up dataframe with sessions but no MNAs yet
    MNAs.w <- MNA.function(data=data, site=sites, web=webs[w], species=species[1], sessions=sessions)
    MNAs.w <- MNAs.w[,1:which(names(MNAs.w)=="Prim")]  
    
    for(i in 1:length(species)){ 
      MNA.sp <- MNA.function(data=data, site=sites, web=webs[w], species=species[i], sessions=sessions)[,ifelse(interpolate==TRUE,7,6)]
      MNAs.w <- data.frame(MNAs.w,MNA.sp)
      names(MNAs.w)[which(names(MNAs.w)=="MNA.sp")] <- species[i]
      
    }
    MNAs.w <- data.frame(site=sites, web=as.character(webs[w]), MNAs.w)
    if(w==1){
      MNAs <- MNAs.w
    }else{
      MNAs <- rbind(MNAs.w,MNAs) 
    }  
    w <- w+1
  }
  # dataframe of just MNAs - removing sessions etc
  mna <- MNAs[,(which(names(MNAs)=="Prim")+1):dim(MNAs)[2]]
  # remove pm from the diversity calculations if indicated
  if(include.pm==FALSE){
    mna <- mna[,-which(names(mna)=="pm")]
  }
  ShannonH <- diversity(mna,index="shannon")
  SimpsonD <- diversity(mna,index="simpson")
  invSimpsonD <- diversity(mna,index="invsimpson")
  speciesN <- specnumber(mna)
  
  # sum up all peromyscus species (other than pm)
  pero.ind <- match(c("pb","pe","pl","pn","ps","pt"),names(MNAs))
  pero.ind <- pero.ind[which(is.finite(pero.ind))]
  peros <- apply(MNAs[,pero.ind],1,sum)
  # sum up all species density besides pm 
  other.sp <- apply(mna,1,sum)
  
  
  MNAs.diversity <- data.frame(MNAs ,ShannonH,SimpsonD,invSimpsonD,speciesN,peros,other.sp)
  
  if(scale==TRUE){
    # leave columns 1-7, the rest of the columns divide by max of column
    scale.MNAs <- apply(MNAs.diversity[,8:dim(MNAs.diversity)[2]],2,function(x){x/max(x,na.rm=TRUE)}) 
    MNAs.diversity <- cbind(MNAs.diversity[,1:7],scale.MNAs)
  }
  
  return(MNAs.diversity)
}