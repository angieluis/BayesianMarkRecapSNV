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
      
     
     # data frame of species MNAs and diversities from
     # diversity.df.function(), species 2 letter codes in lowercase
     # cov.list below must match column headers
     # input should be scaled (scale=TRUE in fxn)
     diversity.data = NULL, 
     
     #  e.g. c("Zuni.1","Zuni.2", "Navajo.1","Navajo.2"
     site.webs = c("Zuni.1", "Zuni.2", "Navajo.1", "Navajo.2",
                   "Grandcanyon.E", "Grandcanyon.M", "Grandcanyon.T"),
   
     # list of temporal covariates and their time lags,
     # e.g, list(ndvi=0,ndvi=1,tmax=3) means use ndvi with no lag
     # and with a lag 1 and tmax with lag 3
     cov.list = list(ndvi = 8, prcp = 0, tmax = 0)
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
     y <- last.year - first.year
     ms <- c(first.month:12, rep(1:12, y - 1), 1:last.month)
     ys <- c(rep(first.year, 
                 length(first.month:12)), 
             rep((first.year + 1):(last.year - 1), 
                 each = 12), 
             rep(last.year, 
                 length(1:last.month)))
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
     
##### ISSUE - datas and data.w not binding
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
       
       # now paste in diversity data if present
       if (length(diversity.data) > 0) {
         # use dates to line up diversity data to temporal data, 
         # then add NAs elswhere
         diversity.data$date <- lubridate::dmy(paste("1", 
                                                     diversity.data$month, 
                                                     diversity.data$year, 
                                                     sep = "-"))
         diversity.data$site <- str_to_lower(diversity.data$site)
         diversity.data$site.web <- paste(diversity.data$site, 
                                          tolower(diversity.data$web))
         datas <- dplyr::left_join(datas, diversity.data)
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
       
       # these columns have NAs that need to be filled for CJS models to run, 
       # they are likely the first time points when there is a lag, so fill in 
       # the next value (this is after removing the first 7 column because that 
       # is just the session numbers etc. not the covariate data)
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
  # }
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