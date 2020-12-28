#### Angie added functions below

## Functions and packages for analysis -------------------------------------- ##


 # packages
   library(lubridate)
   library(daymetr)
   library(reshape2)
   library(zoo)
   library(fitdistrplus)
   library(corrplot)
   library(Hmisc)
   library(PerformanceAnalytics)
   library(R2jags)
   library(mcmcplots)
   library(tidyverse) # dplyr and MASS don't play well (i.e., select and rename)
   library(googledrive)
   
   
   
## Associate GoogleDrive with project --------------------------------------- ##
   
   
 # use drive to store BA RData files that can't be pushed to github
   #drive_auth(email = "ecweid@gmail.com")

 # links to folders 
   # Y8fj7RCdpkgq8QPm67onBCTzARUlmbmB
   # 15rlsr2h9MC9uxOvrSjK1vqD-3eb7UNcj
 # specify file id for ThesisResearch stored on google drive
   #ThesisResearchDrive = as_dribble(as_id("11RiXw-J_1K5W5w5Vmho-ns3LJB8cPkNZ"))
  
  
 # share with angie
   #drive_share(file         = ThesisResearchDrive, 
    #           role         = "writer", 
     #          type         = "user", 
      #         emailAddress = "angela.luis@mso.umt.edu", 
       #        emailMessage = "Let me know if you can access the RData files.")
   
   
 # store the url
   #folder_url <- "https://drive.google.com/drive/folders/11RiXw-J_1K5W5w5Vmho-ns3LJB8cPkNZ"
   
   
 # identify the folder on Drive
   #folder <- drive_get(as_id(folder_url))
   
   
 # identify the files in that folder
   #RData_files <- drive_ls(folder)
   
   
 # download
 # these are saved directly to the ThesisResearch on your computer
   #walk(RData_files$id, ~ drive_download(as_id(.x)))
   
   
   
## Function to check NAs ---------------------------------------------------- ##
   
   
 # function to check, sum, and display all columns with na data
   na.fun <- function(data, ...) {
     data %>%
       dplyr::select(everything()) %>%
       summarise_all(funs(sum(is.na(.))))
   }
   
   
## Function to look for duplicates ------------------------------------------ ##
   
   
   # does not show first occurance
   dupes.fun <- function(data, ...) {
     data %>%
       group_by(.dots = lazyeval::lazy_dots(...)) %>%
       mutate(n = n()) %>%
       filter(n() > 1)
   }
   
   
## Logit function ----------------------------------------------------------- ##
   
   
   logit <- function(x) {
     log(x / (1 - x))
   }
   
   
# Reverse logit funciton ---------------------------------------------------- ##
   
   
   rev.logit <- function(x) {
     exp(x) / (1 + exp(x))
     }

   
## CJS initial values ------------------------------------------------------- ##
   
   
 # function to create matrix of initial values for latent state z
 # we shouldn't give initial values for those elements of z whose value 
 # is specified in the data they get an NA
   cjs.init.z <- function(ch, f) {
     for (i in 1:dim(ch)[1]) {
       if (sum(ch[i, ]) == 1) next
       n2 <- max(which(ch[i, ] == 1))
       ch[i, f[i]:n2] <- NA
       }
     for (i in 1:dim(ch)[1]) {
       ch[i, 1:f[i]] <- NA
       }
     return(ch)
   }
   
   
## Known state values ------------------------------------------------------- ##   
   
   
 # function to create matrix with info about known latent state z
  known.state.cjs <- function(ch) {
    state <- ch
    for (i in 1:dim(ch)[1]) {
      n1 <- min(which(ch[i, ] == 1))
      n2 <- max(which(ch[i, ] == 1))
      state[i, n1:n2] <- 1
      # only filling in those that were 0s but we know were alive
      # because caught before and after
      state[i, n1] <- NA
      }
    state[state == 0] <- NA
    return(state)
  }
  

## Primary CH from secondary CH --------------------------------------------- ##
  
 # function to create primary CH from secondary list
   primary.ch.fun <- function(CH.secondary) { # as list of monthly matrices
     CH.primary <- matrix(NA, 
                          nrow = dim(CH.secondary[[1]])[1], 
                          ncol = length(CH.secondary))
     for (i in 1:dim(CH.primary)[2]) {
       CH.primary[, i] <- apply(CH.secondary[[i]], 1, function(x) {
         ifelse(length(which(is.na(x))) == length(x), NA, sum(x, na.rm = TRUE))
       })
       }
     CH.primary <- replace(CH.primary, CH.primary > 1, 1)
     return(CH.primary)
   }
   
   
## Monthly primary CH ------------------------------------------------------- ##
   
   
 # create monthly capture history (including months not trapped) to pass to the
 # initial values and known state functions make months not trapped = 0 
 # (so not really capture history)
   monthly.primary.CH.fun <- function(CH.primary, temporal.covariates) {
     CH <- matrix(NA, 
                  ncol = max(temporal.covariates$long.month), 
                  nrow = dim(CH.primary)[1])
     old <- which(is.finite(temporal.covariates$Prim))
     CH[, old] <- CH.primary
     CH <- replace(CH, is.na(CH), 0)
     return(CH)
   }
   
   
## Monthly long data from CH ------------------------------------------------ ##
   
   
 # function to take Robust Design capture history list of primary 
 # occasions (months) and turn it into long data frame where each 
 # row is a month (including months not trapped)
  monthly.longdata.CH.fun <- function(CH.secondary,
                                     temporal.covariate.df = NULL,
                                     individual.covariates = NULL,
                                     p_or_c = FALSE) {
    # add the primary occassion to the data (if not already there)
    if (length(which(colnames(CH.secondary[[1]]) == "Prim")) == 0) {
      for (i in 1:length(CH.secondary)) {
        CH.secondary[[i]] <- cbind(data.frame(Prim = i), CH.secondary[[i]])
        }
    }
    
    # define first capture (f) if not already there
    # z will include months not trapped, so first.caught needs to reflect 
    # long.month not primary
    if (length(which(names(individual.covariates) == "f.longmonth")) == 0) {
      CH.primary <- primary.ch.fun(CH.secondary)
      first.caught <- apply(CH.primary, 1, function(x) {
      min(which(x > 0))
        }) # gives primary occasion first caught (not long.month)
      individual.covariates$f.longmonth <- temporal.covariate.df$long.month[match(first.caught, 
                                                                                  temporal.covariate.df$Prim)]
      }
    
    # josh's code to turn into longdata
    obs.dat <- purrr::map_df(
      CH.secondary,
      ~ tibble::as_tibble(.x) %>% 
        dplyr::mutate(
          ID = 1:n()
          ) %>%
        tidyr::gather(Sec, State, -Prim, -ID) %>%
        dplyr::select(ID, Prim, Sec, State)
      )
    
    # join with temporal covariates to make monthly
    obs.dat.full <- inner_join(obs.dat, 
                               temporal.covariate.df[, c("long.month", "Prim")])
    obs.dat.full <- arrange(obs.dat.full, long.month, ID)
    
    #  Subset observation data to observed bits
    obs.dat.full <- dplyr::left_join(obs.dat.full, 
                                     individual.covariates[, c("ID", "f.longmonth")]) %>%
      dplyr::filter(long.month >= f.longmonth)
    return(obs.dat.full)
  }  
  
  
## Multisite p or c array --------------------------------------------------- ##
  
  
 # this accounts for certain webs not being trapped when others were
 # could be used in place of above for not multiple sites
   multisite.p.or.c.array.fun <- function(CH.secondary, # can be list or array
                                          temporal.covariates # df from list
                                          ) {
     Ch.primary <- primary.ch.fun(CH.secondary)
     nind <- dim(CH.secondary[[1]])[1]
     temp.data <- temporal.covariates
     nt <- dim(temp.data)[1]
     n.sec <- unlist(lapply(CH.secondary, function(x) {
       dim(x)[2]
       }))
     p.or.c <- array(NA, dim = c(nind, nt, max(n.sec)))
     
     for (i in 1:nind) {
       
       # sessions this indiv could have been trapped
       m.trapped <- which(is.finite(Ch.primary[i, ]))
       
       # what long months
       longm.trapped <- temp.data$long.month[match(m.trapped, temp.data$Prim)] 
       
       for (m in 1:length(m.trapped)) {
         
         # number of secondary occasions there was trapping just for that 
         # individual(varies by web) that month
         nsec <- length(which(is.finite(CH.secondary[[m.trapped[m]]][i, ]))) 
         for (d in 1:nsec) {
           dsum <- sum(CH.secondary[[m.trapped[m]]][i, 1:(d - 1)])
           p.or.c[i, longm.trapped[m], d] <- ifelse(dsum == 0, 0, 1)
           } # d
         } # m
       } # i
     return(p.or.c)
   }
   
   
## Multisite monthly covariate function ------------------------------------- ##
   

 # funtion to line up covariates to site  
   # updated to include diversity/density data for multiple webs
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
     
     multistate = FALSE, # if multistate model, 
     # need state at first capture
     
     
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
     remove.na = FALSE #
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
     if(length(not.trapped) > 0){
      for (i in 1:length(not.trapped)) {
       session.num[not.trapped[i]] <- session.num[max(which(sn < not.trapped[i]))]
      }
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
         #if (length(fd) == 1) {
         #  data.w <- cbind(data.w, 
         #                  datas[(fd - cov.list[[c]]):(ld - cov.list[[c]]), 
         #                         col])
         #  names(data.w)[dim(data.w)[2]] <- nam
         #}
         #if (length(fd) > 1) { ## if not >1 then prob with '.web' below. I think can do this no matter the length(fd), don't need the bit above
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
         #}
       }
       month.data <- dplyr::left_join(month.data, data.w)
       ind.na <- which(is.na(month.data), arr.ind = TRUE)
       
       # these columns have NAs that need to be filled for CJS models to run, 
       # they are likely the first time points when there is a lag, so fill in 
       # the next value (this is after removing the first 7 column because that 
       # is just the session numbers etc. not the covariate data)
       # May need to consider interpolating instead?
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
     if(multistate==TRUE){
        # for multistate models, need to state at first caught:
        ms.CH.primary <-  primary.MSch.fun(CH.secondary)
        f.state <- numeric()
        for(i in 1:dim(ic)[1]){
          f.state[i] <- ms.CH.primary[i,first.caught[i]]
        }
        ic$f.state <- f.state
     }
     
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
       
       dat <- month.data[, cols] # problem when only 1 column (1 site.web)
       names(dat) <- web.nam
       mat <- matrix(NA, 
                     ncol = dim(month.data)[1], 
                     nrow = dim(individual.covariates)[1])
       for (i in 1:dim(individual.covariates)[1]) {
         w <- individual.covariates$web[i]
         if(length(cols) > 1){ # if more than one site.web
           mat[i, ] <- dat[, which(names(dat) == w)]
         } else{ # if only 1 site.web
           mat[i, ] <- dat
         }
        } #i
       covariate.data[[c + 2]] <- mat
       names(covariate.data)[[c + 2]] <- nam
     } # c
     return(covariate.data)
   }
   
   

   
## Multisite session list function ------------------------------------------ ##
   

 # session list function
  multisite.session.list.fun <- function(
    
    # assuming this data has all the dates 
    # (including those trapped but no animals/pema were caught)
    dirty.data = southwest.dirty, 
    site.webs = c("Zuni.1", "Zuni.2", "Navajo.1", "Navajo.2",
                  "GrandCanyon.E", "GrandCanyon.F", "GrandCanyon.G",
                  "GrandCanyon.M", "GrandCanyon.T")) {
  
    # make everything lowercase
    site.webs <- tolower(site.webs)
    names(dirty.data) <- tolower(names(dirty.data))
    dirty.data$site <- tolower(dirty.data$site)
    dirty.data$web <- tolower(dirty.data$web)
    

    # create site.web column in data
    dirty.data$site.web <- paste(dirty.data$site, dirty.data$web, sep = ".")

    # cut data to just the site.webs wanted
    ind.dirty <- numeric()
    for (i in 1:length(site.webs)) {
      ind.dirty <- c(ind.dirty, which(dirty.data$site.web == site.webs[i]))
      }
    dirty.data <- dirty.data[ind.dirty, ]
    all.sessions <- sort(unique(dirty.data$session))
    session.list <- list(all.sessions = all.sessions)
    for (i in 1:length(site.webs)) {
      session.list[[i + 1]] <- sort(unique(dirty.data$session[which(dirty.data$site.web == site.webs[i])]))
      names(session.list)[i + 1] <- paste("web", site.webs[i], sep = ".")
      }
    return(session.list)
  }
  
  
## Multisite CJS capture history function ----------------------------------- ##
  
  
  multisite.CJS.capture.history.fun <- function(
    
    # assuming this data has all the dates 
    # including those trapped but no animals/pema were caught
    dirty.data = southwest.dirty,
    cleaned.data = southwest.final.clean, # 'cleaned'
    site.webs = c("Zuni.1", "Zuni.2", "Navajo.1", "Navajo.2",
                  "GrandCanyon.E", "GrandCanyon.F", "GrandCanyon.G",
                  "GrandCanyon.M", "GrandCanyon.T"),
    species = "PM") {
    
    # make everything lowercase
    site.webs <- tolower(site.webs)
    species <- tolower(species)
    names(dirty.data) <- tolower(names(dirty.data))
    dirty.data$site <- tolower(dirty.data$site)
    dirty.data$web <- tolower(dirty.data$web)
    names(cleaned.data) <- tolower(names(cleaned.data))
    cleaned.data$site <- tolower(cleaned.data$site)
    cleaned.data$letter_2 <- tolower(cleaned.data$letter_2)
    cleaned.data$web <- tolower(cleaned.data$web)

    # create site.web column in data
    dirty.data$site.web <- paste(dirty.data$site, 
                                 dirty.data$web, 
                                 sep = ".")
    cleaned.data$site.web <- paste(cleaned.data$site, 
                                   cleaned.data$web, 
                                   sep = ".")
    
    # cut data to just the site.webs wanted
    ind.dirty <- numeric()
    ind.clean <- numeric()
    for (i in 1:length(site.webs)) {
      ind.dirty <- c(ind.dirty, which(dirty.data$site.web == site.webs[i]))
      ind.clean <- c(ind.clean, which(cleaned.data$site.web == site.webs[i]))
      }
    dirty.data <- dirty.data[ind.dirty, ]
    cleaned.data <- cleaned.data[ind.clean, ]

    # new date column in right format
    dirty.data$date1 <- as.Date(gsub(" ", 
                                     "", 
                                     dirty.data$date), 
                                format = "%m/%d/%Y")
    cleaned.data$date1 <- as.Date(gsub(" ", 
                                       "", 
                                       cleaned.data$date), 
                                  format = "%m/%d/%Y")

    # new site.tag column of data
    cleaned.data$site.tag <- paste(cleaned.data$site, 
                                   cleaned.data$tag, 
                                   sep = ".")
    all.sessions <- sort(unique(dirty.data$session))
    web.dates <- list()
    
    for (i in 1:length(site.webs)) {
      web.dates[[i]] <- sort(unique(dirty.data$date1[which(dirty.data$site.web == site.webs[i])]))
      }
    names(web.dates) <- paste("web", site.webs, sep = ".")

    # calculate the number of secondary occasions for each sessions at each site
    sec.occ.list <- list()
    for (i in 1:length(site.webs)) {
      
      # unique sessions for that site.web
      s <- sort(unique(dirty.data$session[dirty.data$site.web == site.webs[[i]]])) 
      sec.occ.list[[i]] <- list()
      for (m in 1:length(s)) {
        sec.occ.list[[i]][[m]] <- sort(unique(dirty.data$date1[dirty.data$site.web == site.webs[[i]] & dirty.data$session == s[m]]))
        }
      names(sec.occ.list[[i]]) <- s
      }
    names(sec.occ.list) <- site.webs
    n.sec.occ.list <- list()
    for (i in 1:length(site.webs)) {
      n.sec.occ.list[[i]] <- unlist(lapply(sec.occ.list[[i]], length))
      }
    names(n.sec.occ.list) <- site.webs

    # separate out deermice. These are the site.tag numbers I want to keep (are pema).
    sp.tags <- sort(unique(cleaned.data$site.tag[which(cleaned.data$letter_2 == species)]))
    if (length(which(sp.tags == -9)) > 0) {
      sp.tags <- sp.tags[-which(sp.tags == -9)]
      }
    
    # what are the rows of data that have those tags?
    sp.ind <- numeric()
    for (i in 1:length(sp.tags)) {
      sp.ind <- c(sp.ind, which(cleaned.data$site.tag == sp.tags[i]))
    }
    
    # cleaned data for sites/webs with species in which i'm interested
    sp.data <- cleaned.data[sp.ind, ] 
    
    # max number of secondary occasions per primary (among all sites)
    max.sec.occ <- numeric() 
    for (i in 1:length(all.sessions)) {
      max.sec.occ[i] <- max(unlist(lapply(n.sec.occ.list, function(x) {
      x[which(names(x) == all.sessions[i])]
      })))
    }
    
    # first set up blank capture history list to fill in below
    all.IDs <- sort(unique(sp.data$site.tag))
    Ch.list <- list()
    for (i in 1:length(all.sessions)) {
      mat <- matrix(NA, nrow = length(all.IDs), ncol = max.sec.occ[i])
      rownames(mat) <- all.IDs
      colnames(mat) <- 1:(dim(mat)[2])
      Ch.list[[i]] <- mat
    }
    
    # fill in by web
    for (w in 1:length(site.webs)) {
      web.IDs <- unique(sp.data$site.tag[which(sp.data$site.web == site.webs[w])])
      Session.days <- sec.occ.list[[w]]
      web.dat <- sp.data[which(sp.data$site.web == site.webs[w]), ]
      for (m in 1:length(Session.days)) {
        session <- which(all.sessions == names(Session.days)[m])
        days <- Session.days[[m]]
        for (d in 1:length(days)) {
          for (i in 1:length(web.IDs)) {
            indiv <- which(rownames(Ch.list[[1]]) == web.IDs[i])
            id.dat <- web.dat[which(web.dat$date == days[d] & web.dat$site.tag == web.IDs[i]), ]
            Ch.list[[session]][indiv, d] <- ifelse(dim(id.dat)[1] > 0, 1, 0)
            } # i
          } # d
        cat("web = ", w, "; session = ", m, "\n")
        } # m
      } # w
    return(Ch.list)
  }
  

  
  
  
  
    
## -------------------------------------------------------------------------- ##
### Functions angie added
## -------------------------------------------------------------------------- ##
  
  
## Diversity Functions ------------------------------------------------------ ##  
  
  ################################################################
  # function to calculate MNAs (in dataframe with sessions & 
  # longmonths, including sessions not trapped as NAs)
  ################################################################
  
  library(zoo)
  MNA.function <- function(data=southwest.final.clean, # capture data 
                           site="Zuni", # 1 site at a time
                           web=1, # 1 web at a time
                           species="pm", # can put in more than 1
                           sessions=session.list_Z1$web.zuni.1 # sessions trapped
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
    if(y==0){
      ms <- first.month:last.month
      ys <- rep(first.year,length(ms))
    }
    if(y==1) {
      ms <- c(first.month:12,rep(1:12,y-1),1:last.month)
      ys <- c(rep(first.year,length(first.month:12)),rep(last.year,length(1:last.month)))
    }
    if(y>1){
      ms <- c(first.month:12,rep(1:12,y-1),1:last.month)
      ys <- c(rep(first.year,length(first.month:12)),rep((first.year+1):(last.year-1),each=12),rep(last.year,length(1:last.month)))
    } 
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
          ch[i,j] <- length(which(data$tag==tags[i] & data$session==sessions[j]))
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
    
    
    MNA.data <- data.frame(long.month=1:length(ms),session= s1,year=ys,month=ms,Prim=Prim, MNA=unname(MNA[match(s1,sessions)]), MNA.interp = na.approx(MNA[match(s1,sessions)],1:length(s1)))
    
    return(MNA.data)
    
  }
  
  
  
  
  
  ################################################################
  # function to create dataframe of all species MNAs for a site
  # and diversity indices
  ################################################################
  #Using the vegan package to calculate diversity indices
  library(vegan)
  library(stringr)
  library(dplyr)
  diversity.df.function <- function(
    data=southwest.final.clean, # capture data 
    site="Zuni", # 1 site at a time
    web=1, # 1 web at a time
    sessions=session.list_Z1$web.zuni.1, # 
    interpolate=FALSE, # do you want to interpolate NAs?
    scale=FALSE, # do you want to scale each one between 0 and 1 (divide by max) - this is needed for CJS models (but might want to do it after other data are pasted in)
    include.pm = TRUE # reviewers wanted diversity indices calculated with pm removed (include.pm=FALSE), but removing pm leads to Inf in invSimpsonD calculations (when no other species are present)
    
  ){
    site <- tolower(site)
    web <- tolower(web)
    
    site1 <- site
    web1 <- web
    
    ## moved this up, so all species are included (not just those at this web)
    species <- sort(unique(str_to_lower(data$letter_2)))
    #data <- filter(data, site %in% site, web %in% web )
    data <- filter(data, site == site1, web == web1 )
    
    
    # these are the cleaned data, so no -9 or blank
      #if(length(which(species=="-9"))>0){
      #  species <- species[-which(species=="-9")]
      #}
      #if(length(which(species==""))>0){
      #  species <- species[-which(species=="")]
      #}
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
    #sb <- which(t==0) 
    #if(length(sb)>0){
    #  species <- species[-sb]
    #}   
    
    
    
    w <- 1
    while(w <= length(web)){
      # set up dataframe with sessions but no MNAs yet
      MNAs.w <- MNA.function(data=data, site=site, web=web[w], species=species[1], sessions=sessions)
      MNAs.w <- MNAs.w[,1:which(names(MNAs.w)=="Prim")]  
      
      for(i in 1:length(species)){ # all the species MNAs
        MNA.sp <- MNA.function(data=data, site=site, web=web[w], species=species[i], sessions=sessions)[,ifelse(interpolate==TRUE,7,6)]
        MNAs.w <- data.frame(MNAs.w,MNA.sp)
        names(MNAs.w)[which(names(MNAs.w)=="MNA.sp")] <- species[i]
        
      }
      #### now also calculate MNI
      MNI <- MNA.function(data=data[which(data$snv_pos==1),], site=site, web=web[w], species="pm", sessions=sessions)[,ifelse(interpolate==TRUE,7,6)]
      MNAs.w <- data.frame(MNAs.w,MNI)
      ####
      MNAs.w <- data.frame(site=site, web=as.character(web[w]), MNAs.w)
      if(w==1){
        MNAs <- MNAs.w
      }else{
        MNAs <- rbind(MNAs.w,MNAs) 
      }  
      w <- w+1
    }
    # if only pm in data dimensions are off if do rest. don't need to do it
    if(length(species) > 1){
      # dataframe of just MNAs - removing sessions etc
      mna <- MNAs[,(which(names(MNAs)=="Prim")+1):dim(MNAs)[2]]
      mna <- mna[,-which(names(mna)=="MNI")] #MNI isn't a sep species
    
      # remove pm from the diversity calculations if indicated
      if(include.pm==FALSE){
        mna <- mna[,-which(names(mna)=="pm")]
      }
      ShannonH <- diversity(mna,index="shannon")
      SimpsonD <- diversity(mna,index="simpson")
      invSimpsonD <- diversity(mna,index="invsimpson")
      # if no animals were present, you get Inf, so change to 0
      invSimpsonD <- replace(invSimpsonD,invSimpsonD==Inf,0)
      speciesN <- specnumber(mna)
    
      # sum up all peromyscus species (other than pm)
      pero.ind <- match(c("pb","pe","pl","pn","ps","pt"),names(MNAs))
      pero.ind <- pero.ind[which(is.finite(pero.ind))]
      peros <- apply(MNAs[,pero.ind],1,sum)
    
      # sum up all species density besides pm 
      if(include.pm==TRUE){ # if not removed before, need to remove pm
        mna <- mna[,-which(names(mna)=="pm")]
      }
      othersp <- apply(mna,1,sum) 
    
    
      MNAs.diversity <- data.frame(MNAs ,ShannonH, SimpsonD,
                                   invSimpsonD, speciesN, peros, othersp)
    } else { # if only pm in data, don't need diversity metrics
        MNAs.diversity <- MNAs
    }
    
    if(scale==TRUE){
      # leave columns 1-7, the rest of the columns divide by max of column
      scale.MNAs <- apply(MNAs.diversity[,8:dim(MNAs.diversity)[2]],2,function(x){x/max(x,na.rm=TRUE)}) 
      MNAs.diversity <- cbind(MNAs.diversity[,1:7],scale.MNAs)
    }
    
    return(MNAs.diversity)
  }
  

  
  
  
  
## Functions for Multistate Models --------------------------------------------- ##  
 
  
  
  ############################################################################
  # Robust design Multi-State Capture Histories
  # where states are:
  #   uninfected = "1"
  #   infected = "2"
  #   and optionally those unknown (as -9 in dataset) ="3", otherwise assume negative
  ############################################################################
  
  
  multisite.MS.capture.history.fun <- function(
    # assuming this data has all the dates 
    # including those trapped but no animals/pema were caught
    dirty.data = southwest.dirty,
    cleaned.data = southwest.final.clean, # 'cleaned'
    site.webs = c("Zuni.1", "Zuni.2"),
    species="PM",
    SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
  ){
    
    # make everything lowercase
    site.webs <- tolower(site.webs)
    species <- tolower(species)
    names(dirty.data) <- tolower(names(dirty.data))
    dirty.data$site <- tolower(dirty.data$site)
    dirty.data$web <- tolower(dirty.data$web)
    names(cleaned.data) <- tolower(names(cleaned.data))
    cleaned.data$site <- tolower(cleaned.data$site)
    cleaned.data$letter_2 <- tolower(cleaned.data$letter_2)
    cleaned.data$web <- tolower(cleaned.data$web)
    
    # create site.web column in data
    dirty.data$site.web <- paste(dirty.data$site, 
                                 dirty.data$web, 
                                 sep = ".")
    cleaned.data$site.web <- paste(cleaned.data$site, 
                                   cleaned.data$web, 
                                   sep = ".")
    
    # cut data to just the site.webs wanted
    ind.dirty <- numeric()
    ind.clean <- numeric()
    for (i in 1:length(site.webs)) {
      ind.dirty <- c(ind.dirty, which(dirty.data$site.web == site.webs[i]))
      ind.clean <- c(ind.clean, which(cleaned.data$site.web == site.webs[i]))
    }
    dirty.data <- dirty.data[ind.dirty, ]
    cleaned.data <- cleaned.data[ind.clean, ]
    
    # new date column in right format
    dirty.data$date1 <- as.Date(gsub(" ", 
                                     "", 
                                     dirty.data$date), 
                                format = "%m/%d/%Y")
    cleaned.data$date1 <- as.Date(gsub(" ", 
                                       "", 
                                       cleaned.data$date), 
                                  format = "%m/%d/%Y")
    
    # new site.tag column of data
    cleaned.data$site.tag <- paste(cleaned.data$site, 
                                   cleaned.data$tag, 
                                   sep = ".")
    all.sessions <- sort(unique(dirty.data$session))
    web.dates <- list()
    
    for (i in 1:length(site.webs)) {
      web.dates[[i]] <- sort(unique(dirty.data$date1[which(dirty.data$site.web == site.webs[i])]))
    }
    names(web.dates) <- paste("web", site.webs, sep = ".")
    
    # calculate the number of secondary occasions for each sessions at each site
    sec.occ.list <- list()
    for (i in 1:length(site.webs)) {
      
      # unique sessions for that site.web
      s <- sort(unique(dirty.data$session[dirty.data$site.web == site.webs[[i]]])) 
      sec.occ.list[[i]] <- list()
      for (m in 1:length(s)) {
        sec.occ.list[[i]][[m]] <- sort(unique(dirty.data$date1[dirty.data$site.web == site.webs[[i]] & dirty.data$session == s[m]]))
      }
      names(sec.occ.list[[i]]) <- s
    }
    names(sec.occ.list) <- site.webs
    n.sec.occ.list <- list()
    for (i in 1:length(site.webs)) {
      n.sec.occ.list[[i]] <- unlist(lapply(sec.occ.list[[i]], length))
    }
    names(n.sec.occ.list) <- site.webs
    
    # separate out deermice. These are the site.tag numbers I want to keep (are pema).
    sp.tags <- sort(unique(cleaned.data$site.tag[which(cleaned.data$letter_2 == species)]))
    if (length(which(sp.tags == -9)) > 0) {
      sp.tags <- sp.tags[-which(sp.tags == -9)]
    }
    
    # what are the rows of data that have those tags?
    sp.ind <- numeric()
    for (i in 1:length(sp.tags)) {
      sp.ind <- c(sp.ind, which(cleaned.data$site.tag == sp.tags[i]))
    }
    
    # cleaned data for sites/webs with species in which i'm interested
    sp.data <- cleaned.data[sp.ind, ] 
    
    # max number of secondary occasions per primary (among all sites)
    max.sec.occ <- numeric() 
    for (i in 1:length(all.sessions)) {
      max.sec.occ[i] <- max(unlist(lapply(n.sec.occ.list, function(x) {
        x[which(names(x) == all.sessions[i])]
      })))
    }
    
    # first set up blank capture history list to fill in below
    all.IDs <- sort(unique(sp.data$site.tag))
    Ch.list <- list()
    for (i in 1:length(all.sessions)) {
      mat <- matrix(NA, nrow = length(all.IDs), ncol = max.sec.occ[i])
      rownames(mat) <- all.IDs
      colnames(mat) <- 1:(dim(mat)[2])
      Ch.list[[i]] <- mat
    }
    
    # fill in by web
    for (w in 1:length(site.webs)) {
      web.IDs <- unique(sp.data$site.tag[which(sp.data$site.web == site.webs[w])])
      Session.days <- sec.occ.list[[w]]
      web.dat <- sp.data[which(sp.data$site.web == site.webs[w]), ]
      for (m in 1:length(Session.days)) {
        session <- which(all.sessions == names(Session.days)[m])
        days <- Session.days[[m]]
        for (d in 1:length(days)) {
          for (i in 1:length(web.IDs)) {
            indiv <- which(rownames(Ch.list[[1]]) == web.IDs[i])
            id.day.dat <- web.dat[which(web.dat$date == days[d] & web.dat$site.tag == web.IDs[i]), ]
            id.session.dat <- web.dat[which(web.dat$session == names(Session.days)[m] & web.dat$site.tag == web.IDs[i]), ]
            
            if(dim(id.day.dat)[1] == 0){ # if animal wasn't caught that day
              Ch.list[[session]][indiv, d] <- 0        # put in a 0
            } else{                 # if was caught, check serostatus for all days in that primary period, and put in that status
              status <- id.session.dat$snv_pos
              Ch.list[[session]][indiv, d] <- ifelse(length(which(status==1))>0,2,ifelse(length(which(status==0))>0,1,ifelse(SNV.unknown.state==TRUE,3,1)))
            }
          } #i
        } #d
        cat("web = ", w, "; session = ", m, "\n")
      } #m
    } #w
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

  
  
  ##################################################################
  ## function to create known state z 
  #  1:uninfected
  #  2:infected
  # this works on combined data
  ##################################################################
  known.state.SImsInf <- function(ms=CH.primary.ms){ # ms is multistate primary capture history
    # notseen: label for 'not seen' #here is 3
    state <- ms
    state[state == 0] <- NA
    for(i in 1:dim(ms)[1]){
      n1 <- min(which(ms[i,]>0))
      if(length(which(ms[i, ] == 2)) > 0){ #filling in I's where can
        minI <- min(which(ms[i, ] == 2)) #I's are observation 2
        maxI <- max(which(ms[i, ] == 2))
        state[i, minI:maxI] <- 2}         # I's are state 3
      if(length(which(ms[i, ] == 1)) > 0){  #filling in S's where can
        minS <- min(which(ms[i, ] == 1))  # S's are observation 1
        maxS <- max(which(ms[i, ] == 1))
        state[i, minS:maxS] <- 1}         # S's are state 2
      state[i,n1] <- NA
    }
    
    return(state)
  }
  
  
  
  ##################################################################
  ## function to specify initial values
  # 1 alive as S
  # 2 alive as I
  # 3 dead
  # this just puts in an initial value of alive as last state seen
  ##################################################################
  MSinf.init.z <- function(ch=CH.primary.ms,
                           n.months){ # n.months is the length of months in the dataset by individual because can differ by web [i,m] . I think assumes that all sites start at the same time?
    kn.state <- known.state.SImsInf(ms = ch)
    f <- apply(ch,1,function(x){min(which(x > 0))})
    state <- matrix(NA, nrow = dim(ch)[1], ncol = dim(ch)[2]) 
    # fill in with first state caught
    for(i in 1:dim(ch)[1]){
      f.state <- ch[i,f[i]]
      state[i,] <- rep(f.state,dim(ch)[2]) 
    }
    # remove those that are in the known state
    state <- replace(state,!is.na(kn.state),NA)
    
    for(i in 1:(dim(state)[1])){
      state[i,1:f[i]] <- NA # put NA for when first caught (in likelihood)
      
      if(length(which(kn.state[i,] == 2)) > 0){ 
        maxI <- max(which(kn.state[i,] == 2))
        if(maxI < dim(state)[2] ){
          state[i, (maxI + 1):dim(state)[2]] <- 2 # all after caught as I are I (2)
        }
      }
      if(n.months[i]!=max(n.months)){
        state[i,(n.months[i]+1):dim(ch)[2]] <- NA # replace all after last month in the dataset with NA
      }
    }
    return(state)
  }
  
  
  
##-----------------------------------------------------------------------##

# Functions to combine data from multiple sites 
  
# This assumes that files were created separately for each site (but
  # can contain multiple webs, e.g. an obs.dat for Grand Canyon webs E,M,T,
  # and an obs.dat for Zuni webs 1,2), and we want to paste them together
  # to run several sites in the same model
# IMPORTANT: input data in same order for all functions
# (all these will work for multistate models as well as CJS)
#############################################################
  
  library(rlist)
  
  
  
  #############################################################
  # function to paste all the observation data together (in alphabetical 
  # order) & change IDs to sequential numbering so match with individual 
  # covariates below
  #############################################################
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
  
  
  
  
  
  ######################################################################
  # function to combine all covariates from separate covariate.data lists
  ######################################################################
  
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
  
  
  
  ######################################################################
  # function to combine monthly CHs 
  ######################################################################
  
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
  
  
  ######################################################################
  # function to combine p.or.c matrices 
  ######################################################################
  
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
  
  
  ######################################################################
  # function to combine n.sec.occ 
  ######################################################################

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
  
  
  ######################################################################
  # function to combine months.trapped 
  ######################################################################
  
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

    
  ######################################################################
  # function for cjs.init.z for combined data (because more NAs)
  ######################################################################
  
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
  
  
 