################################################################################
########## Basic Robust Design CJS capture histories (0 or 1) AS LIST
# A matrix for each month (primary occasion), where column is day (secondary occasion) and row is individual. Each matrix has the same number of rows - all individuals have a row each month even though for most they will be zeros.
################################################################################

CJS.capture.history.function <- function(
  data=UNMcaptures, #assuming this data has all the dates (including those trapped but no animals/pema were caught) and is 'cleaned'
  site="Zuni",
  webs=c("1","2"),
  species="PM"
){
  
data$date <- as.Date(gsub(" ", "", data$Date),format="%m/%d/%Y")
dates <- sort(unique(data$date[which(data$site==site)]))
sessions <- sort(unique(data$Session[which(data$site==site)]))

web.dates <- list()
for(i in 1:length(webs)){
  web.dates[[i]] <-sort(unique(data$date[which(data$site==site&data$web==webs[i])]))
}
names(web.dates) <- paste("web",webs,sep=".")

dates.df <- data.frame(dates)#web1=rep("yes",length(dates)),web2=rep("yes",length(Zuni.dates)))
for(i in 1:length(webs)){
  x <- rep("yes", dim(dates.df)[1])
  x[which(is.na(match(dates,web.dates[[i]])))] <- "no"  
  dates.df[,i+1] <- x
  names(dates.df)[i+1] <- paste("web",webs[i],sep=".")
}


### since I want this as a different output maybe make its own function
session.list <- list(all.sessions=sessions) 
for(i in 1:length(webs)){
  session.list[[i+1]] <- sort(unique(data$Session[which(data$site==site&data$web==webs[i])]))
  names(session.list)[i+1] <- paste("web",webs[i],sep=".")
}


time.int <- diff(dates)
first.dates <- dates[c(1,1+which(time.int>1))]

primary.time.int.weeks <- diff(first.dates)/7 # don't think I'm using this now

# data for just these site and webs
inds=numeric()
for(i in 1:length(webs)){
  inds <- c(inds,which(data$site==site & data$web==webs[i]))
}
site.data <- data[inds,] 

# separate out deermice. Since some of the repeaters might have -9 under letter_2, I will go by tag. These are the tag numbers I want to keep (are pema).
sp.tags <- sort(unique(site.data$tag[which(site.data$letter_2==species)]))
sp.tags <- sp.tags[-which(sp.tags==-9)]

# what are the rows of data that have those tags?
sp.ind <- numeric()
for(i in 1:length(sp.tags)){
  sp.ind <- c(sp.ind,which(site.data$tag==sp.tags[i]))
}
site.sp.data <- site.data[sp.ind,] #site&webs I'm interested in, with just the species I'm interested in



Session.days <- list()
for(i in 1:length(sessions)){
  Session.days[[i]] <- sort(unique(site.data$date[which(site.data$Session==sessions[i])]))
    names(Session.days)[i] <- sessions[i]
}


IDs <- sort (unique(site.sp.data$tag))
ID.web <- character()
for(i in 1:length(IDs)){
  ID.web[i] <- as.character(site.sp.data$web[which(site.sp.data$tag==IDs[i])][1])
}


Ch.list <- list()

for(m in 1:length(Session.days)){
  days <- Session.days[[m]]
  ch.mat <- matrix(NA,ncol=length(days),nrow=length(IDs))
  
  for(d in 1:length(days)){
    for(i in 1:length(IDs)){
      col <- which(names(dates.df)==paste("web",ID.web[i],sep="."))
      if(dates.df[which(dates.df$dates==days[d]),col]=="yes"){ # check which web this animal is from and if that web was trapped that day, if not, it will remain an NA
      ch.mat[i,d] <- ifelse(length(which(site.sp.data$tag==IDs[i] & site.sp.data$date==days[d]))>0,1,0)
      }
    }
  }
  dimnames(ch.mat) <- list(IDs,NULL)
  Ch.list[[m]] <- ch.mat
  cat("session = ", m, "\n")
}

return(Ch.list)
 
} 

#####################################################################
############### separate function to return the session list
# list of length # webs + 1
# first is vector of all the sessions for that site
# then a vector for each web, saying which sessions that web
# was trapped
#####################################################################

session.list.function <- function(
  data=UNMcaptures, #assuming this data has all the dates (including those trapped but no animals/pema were caught) and is 'cleaned' and has a column called 'session' or 'Session'
  site="Zuni",
  webs=c("1","2")
){
  
  site=tolower(site)
  names(data) <- tolower(names(data))
  data$site <- tolower(data$site)
  sessions <- sort(unique(data$session[which(data$site==site)]))
  session.list <- list(all.sessions=sessions) 
  for(i in 1:length(webs)){
    session.list[[i+1]] <- sort(unique(data$session[which(data$site==site&data$web==webs[i])]))
    names(session.list)[i+1] <- paste("web",webs[i],sep=".")
  }
  
  return(session.list)
}

############################################################################
# function to create primary CH from secondary list
# (same function as in RobustCJSfunctions.R)
############################################################################
primary.ch.fun <- function(CH.secondary){ # as list of monthly matrices 
  x <- lapply(CH.secondary,rowSums)
  v1 <- unlist(x)
  CH.primary <- matrix(replace(v1, v1>1, 1), nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary)) 
  
  return(CH.primary)
}

############################################################################
# Robust design Multi-State Capture Histories
# where states are:
#   uninfected = "1"
#   infected = "2"
#   and optionally those unknown (as -9 in dataset) ="3", otherwise assume negative
############################################################################

MS.capture.history.function <- function(
  data=UNMcaptures, #assuming this data has all the dates (including those trapped but no animals/pema were caught) and is 'cleaned'
  site="Zuni",
  webs=c("1","2"),
  species="PM",
  SNV.unknown.state=TRUE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
){
  
  data$date <- as.Date(gsub(" ", "", data$Date),format="%m/%d/%Y")
  dates <- sort(unique(data$date[which(data$site==site)]))
  sessions <- sort(unique(data$Session[which(data$site==site)]))
  
  web.dates <- list()
  for(i in 1:length(webs)){
    web.dates[[i]] <-sort(unique(data$date[which(data$site==site&data$web==webs[i])]))
  }
  names(web.dates) <- paste("web",webs,sep=".")
  
  dates.df <- data.frame(dates)#web1=rep("yes",length(dates)),web2=rep("yes",length(Zuni.dates)))
  for(i in 1:length(webs)){
    x <- rep("yes", dim(dates.df)[1])
    x[which(is.na(match(dates,web.dates[[i]])))] <- "no"  
    dates.df[,i+1] <- x
    names(dates.df)[i+1] <- paste("web",webs[i],sep=".")
  }
  
  
  ### since I want this as a different output maybe make its own function
  session.list <- list(all.sessions=sessions) 
  for(i in 1:length(webs)){
    session.list[[i+1]] <- sort(unique(data$Session[which(data$site==site&data$web==webs[i])]))
    names(session.list)[i+1] <- paste("web",webs[i],sep=".")
  }
  
  
  time.int <- diff(dates)
  first.dates <- dates[c(1,1+which(time.int>1))]
  
  primary.time.int.weeks <- diff(first.dates)/7 # don't think I'm using this now
  
  # data for just these site and webs
  inds=numeric()
  for(i in 1:length(webs)){
    inds <- c(inds,which(data$site==site & data$web==webs[i]))
  }
  site.data <- data[inds,] 
  
  # separate out deermice. Since some of the repeaters might have -9 under letter_2, I will go by tag. These are the tag numbers I want to keep (are pema).
  sp.tags <- sort(unique(site.data$tag[which(site.data$letter_2==species)]))
  sp.tags <- sp.tags[-which(sp.tags==-9)]
  
  # what are the rows of data that have those tags?
  sp.ind <- numeric()
  for(i in 1:length(sp.tags)){
    sp.ind <- c(sp.ind,which(site.data$tag==sp.tags[i]))
  }
  site.sp.data <- site.data[sp.ind,] #site&webs I'm interested in, with just the species I'm interested in
  
  
  
  Session.days <- list()
  for(i in 1:length(sessions)){
    Session.days[[i]] <- sort(unique(site.data$date[which(site.data$Session==sessions[i])]))
    names(Session.days)[i] <- sessions[i]
  }
  
  IDs <- sort (unique(site.sp.data$tag))
  ID.web <- character()
  for(i in 1:length(IDs)){
    ID.web[i] <- as.character(site.sp.data$web[which(site.sp.data$tag==IDs[i])][1])
  }
  
  
  Ch.list <- list()
  
  for(m in 1:length(Session.days)){
    days <- Session.days[[m]]
    ch.mat <- matrix(NA,ncol=length(days),nrow=length(IDs))
    
    for(d in 1:length(days)){
      for(i in 1:length(IDs)){
        col <- which(names(dates.df)==paste("web",ID.web[i],sep="."))
        if(dates.df[which(dates.df$dates==days[d]),col]=="yes"){ # check which web this animal is from and if that web was trapped that day, if not, it will remain an NA
          caught <-  ifelse(length(which(site.sp.data$tag==IDs[i] & site.sp.data$date==days[d]))>0,1,0) # caught =1, not caught =0
          if(caught==1){ # if caught, check serostatus for all days in that primary period, and put in that status
            status <- site.sp.data$snv_pos[which(site.sp.data$tag==IDs[i] & site.sp.data$Session==sessions[m])]
            ch.mat[i,d] <- ifelse(length(which(status==1))>0,2,ifelse(length(which(status==0))>0,1,ifelse(SNV.unknown.state==TRUE,3,1)))
          } else{
            ch.mat[i,d] <- 0
          }
          
        }
      }
    }
    dimnames(ch.mat) <- list(IDs,NULL)
    Ch.list[[m]] <- ch.mat
    cat("session = ", m, "\n")
  }
  
  return(Ch.list)
  
} 

####################################################################
## create a primary multi-state capture history
# from a secondary history (collapse to just primary occasions)
####################################################################

primary.MSch.fun <- function(CH.secondary){ # as list of monthly matrices 
  CH.primary <- matrix(NA,nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary))
  for(m in 1:length(CH.secondary)){
    for(i in 1:dim(CH.secondary[[1]])[1]){
      chs <- CH.secondary[[m]][i,]
      if(length(which(is.na(chs)))==length(chs)){
        CH.primary[i,m] <- NA
      } else{
        CH.primary[i,m] <- ifelse(sum(chs)==0,0,unique(chs[which(chs>0)]))
      } 
    }
  }                     
  
  return(CH.primary)
}
