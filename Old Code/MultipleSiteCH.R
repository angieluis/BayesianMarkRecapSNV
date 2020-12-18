### working on robust design capture histories for multiple sites
# working by session with NAs when sessions not trapped at that site

load("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/PemaData.RData")
# this is all data before and after Emily's Cleaning code

# should work with Emily's cleaned and dirty data
dirty.data = pema.dirty # all species all sites before cleaning (not just pema)
cleaned.data = pema.final.clean # all sites/species after cleaning


multisite.CJS.capture.history.function <- function(
  dirty.data = pema.dirty, #assuming this data has all the dates (including those trapped but no animals/pema were caught) 
  cleaned.data = pema.final.clean, # 'cleaned'
  # 2 date columns in these data
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

# need to make sure individual covariates match up now