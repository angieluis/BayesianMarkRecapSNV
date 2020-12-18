
##### recoding covariate data to put in diversity data


Ch.list_GCemt <-  multisite.CJS.capture.history.fun(site.webs=c("grandcanyon.e","grandcanyon.m","grandcanyon.t"))
CH.secondary <- Ch.list_GCemt
tags <- rownames(CH.secondary[[1]])

session.list_GCemt <- multisite.session.list.fun(site.webs=c("grandcanyon.e","grandcanyon.m","grandcanyon.t"))
sessions <- session.list_GCemt$all.sessions

site.webs <- c("grandcanyon.e","grandcanyon.m","grandcanyon.t")

cleaned.data = southwest.final.clean
by.sitetag = TRUE
diversity.data = scaled.MNAs.diversity.longdata
remove.na=FALSE
cov.list = list(ndvi3 = 0, prcp = 0, tmax = 0,pm=0,pm=1,MNI=0,MNI=1,invSimpsonD=0)



cleaned.data = southwest.final.clean # cleaned data
CH.secondary = ms.CH.secondary_GC          # as monthly list
# in CH.secondary row names are tags, tag names line up to CH.secondary
tags = rownames(ms.CH.secondary_GC[[1]])
by.sitetag = TRUE
sessions = session.list_GC$all.sessions
temporal.data = sw.temp.data
diversity.data = scaled.MNAs.diversity.longdata
site.webs = c("grandcanyon.E","grandcanyon.F","grandcanyon.G","grandcanyon.M","grandcanyon.T")
# list of temporal covariates and their time lags
cov.list = list(prcp = 0, prcp3 = 0, prcp6 = 0, prcp12 = 0, 
                ndvi = 0, ndvi3 = 0, ndvi6 = 0, ndvi12 = 0, 
                temp = 0, temp3 = 0, temp6 = 0, temp12 = 0, 
                tmin = 0, tmin3 = 0, tmin6 = 0, tmin12 = 0, 
                tmax = 0, tmax3 = 0, tmax6 = 0, tmax12 = 0, 
                swe  = 0, swe3  = 0, swe6  = 0, swe12  = 0,
                swewinter = 0,
                pm = 0, MNI = 0, invSimpsonD = 0, speciesN = 0,
                peros = 0, othersp = 0, pt = 0) 




#### they way we've been using this is using the lags already calculated: ndvi3=0 instead of ndvi=3


# funtion to line up covariates to site  
multisite.monthly.covariate.fun <- function(
  cleaned.data = southwest.final.clean,
  
  # as monthly list
  CH.secondary = Ch.list,
  
  # tag names that line up to CH.secondary
  tags = rownames(Ch.list[[1]]), 
  
  # if TRUE, then tags are specified by site e.g., "Zuni.1101"
  by.sitetag = TRUE, 
  
  # all sessions to include e.g. 199806 (sessions trapped even no pm caught
  sessions = session.list$all.sessions,
  
  # data frame of monthly temporal data, with either a column
  # called date or yearmon (must include months not trapped),
  # must be in long format like sw.temp.data and have all the
  # sites and all temporal data with no time lags
  temporal.data = sw.temp.data, 
  
  
  # long data frame of species MNAs and diversities 
  # species 2 letter codes in lowercase
  # cov.list below must match column headers
  # input should be scaled (max of any variable=1 across all sites)
  diversity.data = scaled.MNAs.diversity.longdata, # or = NULL
  
  #  e.g. c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"
  site.webs = c("Zuni.1", "Zuni.2", "Navajo.1", "Navajo.2",
                "Grandcanyon.E", "Grandcanyon.M", "Grandcanyon.T"),
  
  # list of temporal covariates and their time lags,
  # e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag
  # and with a lag 1 and tmax with lag 3
  cov.list = list(ndvi = 8, prcp = 0, tmax = 0),
  remove.na = FALSE #,
  #div.cov.list = NULL # #additional diversity covariates specified in same way
) {
  
  names(cleaned.data) <- tolower(names(cleaned.data))
  tags <- tolower(tags)
  site.webs <- tolower(site.webs)
  cleaned.data$site.tag <- paste(tolower(cleaned.data$site), 
                                 tolower(cleaned.data$tag), sep = ".")
  cleaned.data$site.web <- paste(tolower(cleaned.data$site), 
                                 tolower(cleaned.data$web), sep = ".")
  temporal.data$site.web <- paste(tolower(temporal.data$site), 
                                  tolower(temporal.data$web), sep = ".")
  nind <- length(tags)
  ind.clean <- numeric()
  for (i in 1:length(site.webs)) {
    ind.clean <- c(ind.clean, which(cleaned.data$site.web == site.webs[i]))
  }
  cleaned.data <- cleaned.data[ind.clean, ]
  ic <- data.frame(ID = 1:nind, tag = tags)
  site.webi <- character()
  sex <- numeric() # 1 male, 0 female
  for (i in 1:nind) {
    if (by.sitetag == FALSE) {
      ind <- which(cleaned.data$tag == tags[i])
    } else {
      ind <- which(cleaned.data$site.tag == tags[i])
    }
    x <- cleaned.data[ind, ]
    x <- x[order(x$session), ]
    site.webi[i] <- as.character(x$site.web[1])
    sex[i] <- max(x$sex) # they aren't NAs but -9
  }
  ic$web <- factor(site.webi)
  ic$sex <- replace(sex, sex == -9, 0.5) # split the difference for unknown sexes
  sessions.trapped <- sort(unique(sessions))
  first.session <- sessions.trapped[1]
  last.session <- sessions.trapped[length(sessions.trapped)]
  first.montha <- strsplit(as.character(first.session), 
                           split = character(0))[[1]][5:6]
  first.month <- as.numeric(paste(first.montha[1], 
                                  first.montha[2], 
                                  sep = ""))
  first.yeara <- strsplit(as.character(first.session), 
                          split = character(0))[[1]][1:4]
  first.year <- as.numeric(paste(first.yeara[1], 
                                 first.yeara[2], 
                                 first.yeara[3], 
                                 first.yeara[4], sep = ""))
  last.montha <- strsplit(as.character(last.session), 
                          split = character(0))[[1]][5:6]
  last.month <- as.numeric(paste(last.montha[1], 
                                 last.montha[2], 
                                 sep = ""))
  last.yeara <- strsplit(as.character(last.session), 
                         split = character(0))[[1]][1:4]
  last.year <- as.numeric(paste(last.yeara[1], 
                                last.yeara[2], 
                                last.yeara[3], 
                                last.yeara[4], 
                                sep = ""))
  y <- last.year-first.year
  if(y==0){
    ms <- first.month:last.month
    ys <- rep(first.year,length(ms))
  }
  if(y==1) {
    ms <- c(first.month:12,rep(1:12,y-1),1:last.month)
    ys <- c(rep(first.year,length(first.month:12)),
            rep(last.year,length(1:last.month)))
  }
  if(y>1){
    ms <- c(first.month:12,rep(1:12,y-1),1:last.month)
    ys <- c(rep(first.year,length(first.month:12)),
            rep((first.year+1):(last.year-1),each=12),rep(last.year,
                                              length(1:last.month)))
  } 
  # yearmonth <- paste(as.character(ms),as.character(ys),sep="")
  mc <- as.character(ms)
  for (i in 1:length(ms)) {
    mc[i] <- ifelse(ms[i] < 10, paste(paste(as.character(ys[i]), 
                                            "0", 
                                            as.character(ms[i]), 
                                            sep = "")), 
                    paste(paste(as.character(ys[i]), 
                                as.character(ms[i]), 
                                sep = "")))
  }
  s1 <- sessions.trapped[match(mc, sessions.trapped)]
  not.trapped <- which(is.na(s1))
  sn <- 1:length(sessions.trapped)
  session.num <- sn[match(mc, sessions.trapped)]
  for (i in 1:length(not.trapped)) {
    session.num[not.trapped[i]] <- session.num[max(which(sn < not.trapped[i]))]
  }
  Prim <- session.num
  Prim[which(is.na(s1))] <- NA
  month.data <- data.frame(long.month = 1:length(ms), 
                           session = s1, 
                           year = ys, 
                           month = ms, 
                           covariate.prim = session.num, 
                           Prim = Prim)
  month.data$season <- ifelse(month.data$month == 12 | 
                                month.data$month == 1 | 
                                month.data$month == 2 | 
                                month.data$month == 3, 1, 
                              ifelse(month.data$month == 4 | 
                                       month.data$month == 5 | 
                                       month.data$month == 6, 2, 
                                     ifelse(month.data$month == 7 | 
                                              month.data$month == 8 | 
                                              month.data$month == 9, 3, 4)))
  
  
  if (length(temporal.data) > 0) { # |length(diversity.data)>0
    ls <- length(site.webs)
    datas <- temporal.data[which(temporal.data$site.web == site.webs[1]), ]
    if (ls > 1) {
      for (i in 2:ls) {
        datas <- rbind(datas, temporal.data[which(temporal.data$site.web == site.webs[i]), ])
      }
    }
    
    # make a wide data frame with date/yearmon going from first trapped 
    # session to last trapped session
    wdate <- lubridate::dmy(paste("1", 
                                  month.data$month, 
                                  month.data$year, 
                                  sep = "-"))
    data.w <- data.frame(date = sort(unique(wdate)))
    data.w$year <- lubridate::year(data.w$date)
    data.w$month <- lubridate::month(data.w$date)
    datas$date <- lubridate::dmy(paste("1", datas$date))

    ##### Now paste in diversity data if present
    if (length(diversity.data) > 0) {
      
      diversity.data$date <- lubridate::ymd(paste(diversity.data$year,diversity.data$month,"1"))
      diversity.data$site.web <- paste(diversity.data$site,diversity.data$web,sep=".")
      diversity.data$site <- as.character(diversity.data$site)
      diversity.data$web <- as.character(diversity.data$web)

      # reduce diversity data to just these site.webs
      ls <- length(site.webs)
      div.data <- diversity.data[which(diversity.data$site.web == site.webs[1]), ]
      if (ls > 1) {
        for (i in 2:ls) {
          div.data <- rbind(div.data, diversity.data[which(diversity.data$site.web == site.webs[i]), ])
        }
      }
        

      datas <- dplyr::left_join(datas, div.data) 
    }
    
                                
    cl <- length(cov.list)
    for (c in 1:cl) {
      nam <- paste(names(cov.list)[c], cov.list[[c]], sep = "_")
      col <- which(names(datas) == names(cov.list)[c])
      fd <- which(datas$date == data.w$date[1])
      ld <- which(datas$date == data.w$date[length(data.w$date)])
      if (length(fd) == 1) {
        data.w <- cbind(data.w, 
                        datas[(fd - cov.list[[c]]):(ld - cov.list[[c]]), 
                              col])
        names(data.w)[dim(data.w)[2]] <- nam
      }
      if (length(fd) > 1) {
        for (i in 1:length(fd)) {
          data.w <- cbind(data.w, 
                          datas[(fd[i] - cov.list[[c]]):(ld[i] - cov.list[[c]]), 
                                col])
          names(data.w)[dim(data.w)[2]] <- paste(nam, 
                                                 paste("web", 
                                                       datas$site.web[fd[i]], 
                                                       sep = ""), 
                                                 sep = ".")
        }
      }
    }
    month.data <- dplyr::left_join(month.data, data.w)
    ind.na <- which(is.na(month.data), arr.ind = TRUE)
    
    #----------------------------------------------------------- What's going on here?
                              # I don't think I need to remove NAs
                              # and shouldn't just put in last value when there's a big gap
    
    # these columns have NAs that need to be filled for CJS models to run, 
    # they are likely the first time points when there is a lag, so fill in 
    # the next value (this is after removing the first 7 column because that 
    # is just the session numbers etc. not the covariate data)
    if(remove.na==TRUE){
      ind.na <- ind.na[-which(ind.na[, 2] < 9), ]
      if (dim(ind.na)[1] > 0) {
        for (i in 1:dim(ind.na)[1]) {
          # find the next finite value after this one and plug it in
          vi <- which(is.finite(month.data[, ind.na[i, 2]]))
          vin <- vi[min(which(vi > ind.na[i, 1]))]
          month.data[ind.na[i, 1], ind.na[i, 2]] <- month.data[vin, ind.na[i, 2]]
        }
      }
    }
    #-----------------------------------------------------------
    
  } #if covariate data
  
  ## add first capture to individual covariates data
  CH.primary <- primary.ch.fun(CH.secondary)
  first.caught <- apply(CH.primary, 1, function(x) {
    min(which(x > 0))
  }) # gives primary occasion first caught (not long.month)
  ic$f.longmonth <- month.data$long.month[match(first.caught, month.data$Prim)]
  individual.covariates <- ic
  covariate.data <- list(individual.covariates = individual.covariates, 
                         temporal.covariates = month.data)
  # make matrices of temporal covariates matched up to individuals (temporal 
  # data for the individual based on which web they were on): dimensions 
  # [i,m], where m is longmonth - all months not just those trapped
  for (c in 1:length(cov.list)) {
    nam <- paste(names(cov.list)[c], cov.list[[c]], sep = "_")
    cols <- grep(nam, names(month.data))
    col.name <- names(month.data)[cols]
    web.nam <- unlist(lapply(strsplit(col.name, ".web"), function(x) {
      x[2]
    }))
    dat <- month.data[, cols]
    names(dat) <- web.nam
    mat <- matrix(NA, 
                  ncol = dim(month.data)[1], 
                  nrow = dim(individual.covariates)[1])
    for (i in 1:dim(individual.covariates)[1]) {
      w <- individual.covariates$web[i]
      mat[i, ] <- dat[, which(names(dat) == w)]
    } # i
    covariate.data[[c + 2]] <- mat
    names(covariate.data)[[c + 2]] <- nam
  } # c
  return(covariate.data)
}



##################################################################
## function to specify initial values
# 1 alive as S
# 2 alive as I
# 3 dead
# this just puts in an initial value of alive as last state seen
##################################################################
MSinf.init.z <- function(ch=CH.primary.ms){
  kn.state <- known.state.SImsInf(ms = ch)
  f <- apply(ch,1,function(x){min(which(x > 0),na.rm=TRUE)})
  last.month.trapped <- apply(ch,1,function(x){max(which(is.finite(x)))})
  
  state <- matrix(1, nrow = dim(ch)[1], ncol = dim(ch)[2]) # default is S (1)
  state <- replace(state,!is.na(kn.state),NA)
  
  for(i in 1:(dim(state)[1])){
    state[i,1:f[i]] <- NA
    
    if(length(which(kn.state[i,] == 2)) > 0){
      maxI <- max(which(kn.state[i,] == 2))
      if(maxI < dim(state)[2] ){
        state[i, (maxI + 1):dim(state)[2]] <- 2 # all after caught as I are I (2)
      }
    }
    state[i,(last.month.trapped[i]+1):dim(state)[2]] <- NA
  }
  return(state)
}



#######################################
## problems here ... returning 0 when shouldn't 


primary.MSch.fun <- function(CH.secondary){ # as list of monthly matrices 
  CH.primary <- matrix(NA,nrow=dim(CH.secondary[[1]])[1], ncol=length(CH.secondary))
  for(m in 1:length(CH.secondary)){
    for(i in 1:dim(CH.secondary[[1]])[1]){
      chs <- CH.secondary[[m]][i,]
      if(length(which(is.na(chs))) == length(chs)){
        CH.primary[i,m] <- NA
      } else{
        chs <- chs[is.finite(chs)]
        CH.primary[i,m] <- ifelse(sum(chs) == 0, 0, unique(chs[which(chs > 0)]))
      } 
    }
  }                     
  
  return(CH.primary)
}