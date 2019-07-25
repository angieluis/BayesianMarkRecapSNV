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


  ### also have this as its own function
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
# Robust Design CJS capture histories (like above) for multiple sites
############################################################################


multisite.CJS.capture.history.function <- function(
  dirty.data = pema.dirty, #assuming this data has all the dates (including those trapped but no animals/pema were caught) 
  cleaned.data = pema.final.clean, # 'cleaned'
  site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"),
  species="PM"
){
  
  # make everything lowercase
  site.webs <- tolower(site.webs)
  species <- tolower(species)
  names(dirty.data) <- tolower(names(dirty.data))
  dirty.data$site <- tolower(dirty.data$site)
  names(cleaned.data) <- tolower(names(cleaned.data))
  cleaned.data$site <- tolower(cleaned.data$site)
  cleaned.data$letter_2 <- tolower(cleaned.data$letter_2)
  
  # create site.web column in data
  dirty.data$site.web <- paste(dirty.data$site,dirty.data$web,sep=".")
  cleaned.data$site.web <- paste(cleaned.data$site,cleaned.data$web,sep=".")
  
  # cut data to just the site.webs wanted
  ind.dirty <- numeric()
  ind.clean <- numeric()
  for(i in 1:length(site.webs)){
    ind.dirty <- c(ind.dirty,which(dirty.data$site.web==site.webs[i]))
    ind.clean <- c(ind.clean,which(cleaned.data$site.web==site.webs[i]))
  }
  dirty.data <- dirty.data[ind.dirty,]
  cleaned.data <- cleaned.data[ind.clean,]
  
  # new date column in right format
  dirty.data$date1 <- as.Date(gsub(" ", "", dirty.data$date),format="%m/%d/%Y")
  cleaned.data$date1 <- as.Date(gsub(" ", "", cleaned.data$date),format="%m/%d/%Y")
  
  # new site.tag column of data
  cleaned.data$site.tag <- paste(cleaned.data$site,cleaned.data$tag,sep=".")
  
  all.sessions <- sort(unique(dirty.data$session))
  
  web.dates <- list()
  for(i in 1:length(site.webs)){
    web.dates[[i]] <-sort(unique(dirty.data$date1[which(dirty.data$site.web==site.webs[i])]))
  }
  names(web.dates) <- paste("web",site.webs,sep=".")
  
  # calculate the number of secondary occasions for each sessions at each site
  sec.occ.list <- list()
  for(i in 1:length(site.webs)){
    s <- sort(unique(dirty.data$session[dirty.data$site.web==site.webs[[i]]]))#unique sessions for that site.web
    sec.occ.list[[i]] <- list()
    for(m in 1:length(s)){
      sec.occ.list[[i]][[m]] <- sort(unique(dirty.data$date1[dirty.data$site.web==site.webs[[i]]&dirty.data$session==s[m]]))
    }
    names(sec.occ.list[[i]]) <- s
  }
  names(sec.occ.list) <- site.webs
  
  n.sec.occ.list <- list()
  for(i in 1:length(site.webs)){
    n.sec.occ.list[[i]] <- unlist(lapply(sec.occ.list[[i]],length))
  }
  names(n.sec.occ.list) <- site.webs
  
  
  # separate out deermice. These are the site.tag numbers I want to keep (are pema).
  sp.tags <- sort(unique(cleaned.data$site.tag[which(cleaned.data$letter_2==species)]))
  if(length(which(sp.tags==-9))>0){
    sp.tags <- sp.tags[-which(sp.tags==-9)]
  }
  
  # what are the rows of data that have those tags?
  sp.ind <- numeric()
  for(i in 1:length(sp.tags)){
    sp.ind <- c(sp.ind,which(cleaned.data$site.tag==sp.tags[i]))
  }
  sp.data <- cleaned.data[sp.ind,] # cleaned data for sites&webs I'm interested in, with just the species I'm interested in
  
  max.sec.occ <- numeric() # max number of secondary occasions per primary (among all sites)
  for(i in 1:length(all.sessions)){
    max.sec.occ[i] <- max(unlist(lapply(n.sec.occ.list,function(x){x[which(names(x)==all.sessions[i])]})))
  }
  
  # first set up blank capture history list to fill in below
  all.IDs <- sort(unique(sp.data$site.tag))
  Ch.list <- list()
  for(i in 1:length(all.sessions)){
    mat <- matrix(NA,nrow=length(all.IDs),ncol=max.sec.occ[i])
    rownames(mat) <- all.IDs
    colnames(mat) <- 1:(dim(mat)[2])
    Ch.list[[i]] <- mat
  }
  
  # fill in by web
  for(w in 1:length(site.webs)){
    web.IDs <- unique(sp.data$site.tag[which(sp.data$site.web==site.webs[w])])
    Session.days <- sec.occ.list[[w]]
    web.dat <- sp.data[which(sp.data$site.web==site.webs[w]),]
    for(m in 1:length(Session.days)){
      session <- which(all.sessions==names(Session.days)[m])
      days <- Session.days[[m]]
      for(d in 1:length(days)){
        for(i in 1:length(web.IDs)){
          indiv <- which(rownames(Ch.list[[1]])==web.IDs[i])
          id.dat <- web.dat[which(web.dat$date==days[d] & web.dat$site.tag==web.IDs[i]),]
          Ch.list[[session]][indiv,d] <- ifelse(dim(id.dat)[1]>0, 1, 0)
        } # i
      } # d
      cat("web = ", w, "; session = ", m, "\n")
    } # m 
  } # w
  
  return(Ch.list)
  
} 


multisite.session.list.function <- function(
      dirty.data = pema.dirty, #assuming this data has all the dates (including those trapped but no animals/pema were caught) 
      site.webs=c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2")
){
  # make everything lowercase
  site.webs <- tolower(site.webs)
  names(dirty.data) <- tolower(names(dirty.data))
  dirty.data$site <- tolower(dirty.data$site)

  # create site.web column in data
  dirty.data$site.web <- paste(dirty.data$site,dirty.data$web,sep=".")

  # cut data to just the site.webs wanted
  ind.dirty <- numeric()
  for(i in 1:length(site.webs)){
    ind.dirty <- c(ind.dirty,which(dirty.data$site.web==site.webs[i]))
  }
  dirty.data <- dirty.data[ind.dirty,]

  all.sessions <- sort(unique(dirty.data$session))
  
  session.list <- list(all.sessions=all.sessions) 
  for(i in 1:length(site.webs)){
    session.list[[i+1]] <- sort(unique(dirty.data$session[which(dirty.data$site.web==site.webs[i])]))
    names(session.list)[i+1] <- paste("web",site.webs[i],sep=".")
  }
  
  return(session.list)
  
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
