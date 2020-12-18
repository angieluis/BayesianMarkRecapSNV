## Format capture histories ------------------------------------------------- ##

## code to run all sites and models for variable selection
## UPDATED - 11/19/2019 - commented out the 18-month time lag to reflect its
## removal from the cov.ind model run. See justification in Updated README

## -------------------------------------------------------------------------- ##


 # load capture data (clean and dirty southwest data)
   load("Data/AllCaptureData.RData")


 # source capture history functions and RDCJS functions
   source("Code/01_sw_data_functions.R")


 # read in normalized covariates
   sw.temp.data <- read.csv("Data/Updated/updated_southwest_covariates_norm.csv")
   
   
## Format RData for Zuni ---------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
   
   
 # use dirty for ALL dates trapped, clean for everything else
 # lower or upper case not important for species or site
   CH.secondary <- multisite.CJS.capture.history.fun(
     dirty.data   = southwest.dirty, 
     cleaned.data = southwest.final.clean,
     site.webs    = c("Zuni.1", "Zuni.2"), 
     species      = "PM")
   
    
 # make a list of every session trapped from dirty data
 # list of 3 vectors - length of web plus one
 # the first element of the list is all the sessions for all of the sites
 # all other elements are sessions trapped for that specific web
   session.list <- multisite.session.list.fun(
     dirty.data = southwest.dirty,
     site.webs  = c("Zuni.1", "Zuni.2")
   )
   
   
 # makes a primary capture history from secondary occasions
   CH.primary <- primary.ch.fun(CH.secondary)
   
   
 # rename secondary occasions from 1 to number of days
 # looking at column names and changing them
   for (m in 1:length(CH.secondary)) {
     colnames(CH.secondary[[m]]) <- 1:dim(CH.secondary[[m]])[2]
   } 
   
   
 # number of secondary occasions - maximum of secondary occasions across webs
   n.sec.occ <- unlist(lapply(CH.secondary, function(x) {
     dim(x)[2]
   }))
   
   
 # apply covariate function to normalized data
   covariate.data <- multisite.monthly.covariate.fun(
     cleaned.data = southwest.final.clean, # cleaned data
     CH.secondary = CH.secondary,          # as monthly list
     
   # in CH.secondary row names are tags, tag names line up to CH.secondary
     tags = rownames(CH.secondary[[1]]),
     
   # if TRUE, then tags are specified by site e.g., "Zuni.1101"
     by.sitetag = TRUE, 
     
   # all sessions to include e.g. 199806 (sessions trapped even if no pm caught)
     sessions = session.list$all.sessions,
     
   # data frame of monthly temporal data, with either a column
   # sw.temp.data dataframe made in covariate script
     temporal.data = sw.temp.data, 
     
   # called date or yearmon (must include months not trapped),
   # must be in long format, all sites, all temporal data, no time lags
     site.webs = c("Zuni.1", "Zuni.2"),
     
   # list of temporal covariates and their time lags
   # this means use ndvi for current month, and ndvi rolling mean with windows 
   # of current through t-3, current through t-6, current through t-12, and 
   # current through t-18, (no lags on these because already in dataframe)
     cov.list = list(prcp = 0, prcp3 = 0, prcp6 = 0, prcp12 = 0, #prcp18 = 0,
                     ndvi = 0, ndvi3 = 0, ndvi6 = 0, ndvi12 = 0, #ndvi18 = 0,
                     temp = 0, temp3 = 0, temp6 = 0, temp12 = 0, #temp18 = 0,
                     tmin = 0, tmin3 = 0, tmin6 = 0, tmin12 = 0, #tmin18 = 0,
                     tmax = 0, tmax3 = 0, tmax6 = 0, tmax12 = 0, #tmax18 = 0,
                     swe  = 0, swe3  = 0, swe6  = 0, swe12  = 0,
                     swewinter = 0) #swe18  = 0)
   ) 
   
   
 # output is a list with the first element individual.covariates data frame
 # second element called temporal.covariates to be used for long capture 
 # histories below other elements of the list are matrices of temporal 
 # covariates matching up individuals by web to their temporal data, e.g., 
 # temporal.covariates$ndvi_2 is a matrix of NDVI lag 2 with dimension 
 # [individual, month] (where month is longmonth not just months trapped), 
 # so can be used for phi~ndvi_0[i,m].
   
   
 # array of dim [indiv,total months, max days] #  those not trapped=NA
 # matches up to p[i,m,d]
   p.or.c <- multisite.p.or.c.array.fun(CH.secondary, 
                                        covariate.data$temporal.covariates)
   
   
## Convert to long format --------------------------------------------------- ##
   
   
 # Josh's code that formats long data
   obs.dat <- monthly.longdata.CH.fun(CH.secondary, 
                                      covariate.data$temporal.covariates, 
                                      covariate.data$individual.covariates)
   
   
 # turns secondary into primary and now we add the months that weren't trapped
 # we need this to fill in the known state and initial value functions
   monthlyCH <- monthly.primary.CH.fun(CH.primary,
                                       covariate.data$temporal.covariates)
   
   
 # accounts for how certain webs weren't trapped certain months, that are
 # within the long data set (i.e., shit wasn't trapped at the same time - month)
 # we use all the covariates but we only trapped some of those times, need
 # covariates for the observations process, simulating state process over all
 # months, even those not trapped
   webmonths <- list()
   
   
 # loop through sessions   
   for (i in 1:(length(session.list) - 1)) {
     x <- match(session.list[[i + 1]], 
                covariate.data$temporal.covariates$session)
     webmonths[[i]] <- x[which(is.finite(x))]
   }
   
   names(webmonths) <- names(session.list)[-1]
   
   
 # matrix for months trapped
   months.trapped.mat <- matrix(NA, 
                                nrow = dim(CH.secondary[[1]])[1], 
                                ncol = max(unlist(lapply(webmonths, length))))
   
   
 # create numeric
   length.months.trapped <- numeric()
   
   
 # loop through months trapped
   for (i in 1:dim(months.trapped.mat)[1]) {
     # this is a factor currently
     webnam <- covariate.data$individual.covariates$web[i] 
     webi   <- which(names(webmonths) == paste("web", webnam, sep = "."))
     length.months.trapped[i] <- length(webmonths[[webi]])
     months.trapped.mat[i, 1:length.months.trapped[i]] <- webmonths[[webi]]
   }
   
   
## Arrange covariates in an array ------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
   
 # put all the different ndvi covariates into an array for model selection,
 # then will estimate an indicator which will give prob that covariate should
 # be included in the model the only thing to change is covariate data above 
 # (i.e, sub temp for ndvi)
   
   
 # create array for prcp
   prcp.covariate.array <- array(c(
     covariate.data$prcp_0, # 1
     covariate.data$prcp3_0, # 2
     covariate.data$prcp6_0, # 3
     covariate.data$prcp12_0 # 4
     #covariate.data$prcp18_0 # 5
   ), dim = c(dim(covariate.data$prcp_0), 4)) #5))
   
   
 # create array for ndvi
   ndvi.covariate.array <- array(c(
     covariate.data$ndvi_0, # 1
     covariate.data$ndvi3_0, # 2
     covariate.data$ndvi6_0, # 3
     covariate.data$ndvi12_0 # 4
     #covariate.data$ndvi18_0 # 5
   ), dim = c(dim(covariate.data$ndvi_0), 4)) #5))
   
   
 # create array for temp
   temp.covariate.array <- array(c(
     covariate.data$temp_0, # 1
     covariate.data$temp3_0, # 2
     covariate.data$temp6_0, # 3
     covariate.data$temp12_0 # 4
     #covariate.data$temp18_0 # 5
   ), dim = c(dim(covariate.data$temp_0), 4)) #5))
   
   
 # create array for tmin
   tmin.covariate.array <- array(c(
     covariate.data$tmin_0, # 1
     covariate.data$tmin3_0, # 2
     covariate.data$tmin6_0, # 3
     covariate.data$tmin12_0 # 4
     #covariate.data$tmin18_0 # 5
   ), dim = c(dim(covariate.data$tmin_0), 4)) #5))
   
   
 # create array for tmax
   tmax.covariate.array <- array(c(
     covariate.data$tmax_0, # 1
     covariate.data$tmax3_0, # 2
     covariate.data$tmax6_0, # 3
     covariate.data$tmax12_0 # 4
     #covariate.data$tmax18_0 # 5
   ), dim = c(dim(covariate.data$tmax_0), 4)) #5))
   
   
 # create array for swe
   swe.covariate.array <- array(c(
     covariate.data$swe_0, # 1
     covariate.data$swe3_0, # 2
     covariate.data$swe6_0, # 3
     covariate.data$swe12_0, # 4
     #covariate.data$swe18_0 # 5
     covariate.data$swewinter_0 # 5
   ), dim = c(dim(covariate.data$swe_0), 5))
   
   
   # take care of problem of differing number of secondary occasions
   n.sec.occ <- matrix(NA,
                       nrow = dim(CH.primary)[1],
                       ncol = dim(CH.primary)[2])
   
   for(t in 1:dim(n.sec.occ)[2]) {
     n.sec.occ[,t] <- apply(CH.secondary[[t]],
                            1,
                            function(x){length(which(is.finite(x)))})
   }
   
   
 # save Zuni data
   save.image("Data/Updated/Zuni_CaptureHistories_Jan2020.RData")
   
   
## Format data for Navajo --------------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
 # clear environment and reload necessary files
 # see above for comments
 # 2/11/2019 - updated covariate arrays with swe.winter
   
   
   CH.secondary <- multisite.CJS.capture.history.fun(
     dirty.data   = southwest.dirty, 
     cleaned.data = southwest.final.clean,
     site.webs    = c("Navajo.1", "Navajo.2"), 
     species      = "PM")
   

   session.list <- multisite.session.list.fun(
     dirty.data = southwest.dirty,
     site.webs  = c("Navajo.1", "Navajo.2")
   )
   
   
   CH.primary <- primary.ch.fun(CH.secondary)
   

   for (m in 1:length(CH.secondary)) {
     colnames(CH.secondary[[m]]) <- 1:dim(CH.secondary[[m]])[2]
   }
   

   n.sec.occ <- unlist(lapply(CH.secondary, function(x) {
     dim(x)[2]
   }))
   

   covariate.data <- multisite.monthly.covariate.fun(
     cleaned.data = southwest.final.clean, # cleaned data
     CH.secondary = CH.secondary,          # as monthly list
     tags = rownames(CH.secondary[[1]]),
     by.sitetag = TRUE, 
     sessions = session.list$all.sessions,
     temporal.data = sw.temp.data,
     site.webs = c("Navajo.1", "Navajo.2"),
     cov.list = list(prcp = 0, prcp3 = 0, prcp6 = 0, prcp12 = 0, #prcp18 = 0,
                     ndvi = 0, ndvi3 = 0, ndvi6 = 0, ndvi12 = 0, #ndvi18 = 0,
                     temp = 0, temp3 = 0, temp6 = 0, temp12 = 0, #temp18 = 0,
                     tmin = 0, tmin3 = 0, tmin6 = 0, tmin12 = 0, #tmin18 = 0,
                     tmax = 0, tmax3 = 0, tmax6 = 0, tmax12 = 0, #tmax18 = 0,
                     swe  = 0, swe3  = 0, swe6  = 0, swe12  = 0,
                     swewinter = 0) #swe18  = 0)
   )
   

   p.or.c <- multisite.p.or.c.array.fun(CH.secondary, 
                                        covariate.data$temporal.covariates)
   
   
## Convert to long format --------------------------------------------------- ##
   
   
   obs.dat <- monthly.longdata.CH.fun(CH.secondary, 
                                      covariate.data$temporal.covariates, 
                                      covariate.data$individual.covariates)
   

   monthlyCH <- monthly.primary.CH.fun(CH.primary,
                                       covariate.data$temporal.covariates)
   

   webmonths <- list()
   
 
   for (i in 1:(length(session.list) - 1)) {
     x <- match(session.list[[i + 1]], 
                covariate.data$temporal.covariates$session)
     webmonths[[i]] <- x[which(is.finite(x))]
   }
   
   names(webmonths) <- names(session.list)[-1]
   

   months.trapped.mat <- matrix(NA, 
                                nrow = dim(CH.secondary[[1]])[1], 
                                ncol = max(unlist(lapply(webmonths, length))))
   

   length.months.trapped <- numeric()
   

   for (i in 1:dim(months.trapped.mat)[1]) {
     webnam <- covariate.data$individual.covariates$web[i] 
     webi   <- which(names(webmonths) == paste("web", webnam, sep = "."))
     length.months.trapped[i] <- length(webmonths[[webi]])
     months.trapped.mat[i, 1:length.months.trapped[i]] <- webmonths[[webi]]
   }
   
   
## Arrange covariates in an array ------------------------------------------- ##
## -------------------------------------------------------------------------- ##  

   
 # create array for prcp
   prcp.covariate.array <- array(c(
     covariate.data$prcp_0, # 1
     covariate.data$prcp3_0, # 2
     covariate.data$prcp6_0, # 3
     covariate.data$prcp12_0 # 4
     #covariate.data$prcp18_0 # 5
   ), dim = c(dim(covariate.data$prcp_0), 4)) #5))
   
   
 # create array for ndvi
   ndvi.covariate.array <- array(c(
     covariate.data$ndvi_0, # 1
     covariate.data$ndvi3_0, # 2
     covariate.data$ndvi6_0, # 3
     covariate.data$ndvi12_0 # 4
     #covariate.data$ndvi18_0 # 5
   ), dim = c(dim(covariate.data$ndvi_0), 4)) #5))
   
   
 # create array for temp
   temp.covariate.array <- array(c(
     covariate.data$temp_0, # 1
     covariate.data$temp3_0, # 2
     covariate.data$temp6_0, # 3
     covariate.data$temp12_0 # 4
     #covariate.data$temp18_0 # 5
   ), dim = c(dim(covariate.data$temp_0), 4)) #5))
   
   
 # create array for tmin
   tmin.covariate.array <- array(c(
     covariate.data$tmin_0, # 1
     covariate.data$tmin3_0, # 2
     covariate.data$tmin6_0, # 3
     covariate.data$tmin12_0 # 4
     #covariate.data$tmin18_0 # 5
   ), dim = c(dim(covariate.data$tmin_0), 4)) #5))
   
   
 # create array for tmax
   tmax.covariate.array <- array(c(
     covariate.data$tmax_0, # 1
     covariate.data$tmax3_0, # 2
     covariate.data$tmax6_0, # 3
     covariate.data$tmax12_0 # 4
     #covariate.data$tmax18_0 # 5
   ), dim = c(dim(covariate.data$tmax_0), 4)) #5))
   
   
 # create array for swe
   swe.covariate.array <- array(c(
     covariate.data$swe_0, # 1
     covariate.data$swe3_0, # 2
     covariate.data$swe6_0, # 3
     covariate.data$swe12_0, # 4
     #covariate.data$swe18_0 # 5
     covariate.data$swewinter_0
   ), dim = c(dim(covariate.data$swe_0), 5)) #5))
   
   
 # take care of problem of differing number of secondary occasions
   n.sec.occ <- matrix(NA,
                       nrow = dim(CH.primary)[1],
                       ncol = dim(CH.primary)[2])
   
   for(t in 1:dim(n.sec.occ)[2]) {
     n.sec.occ[,t] <- apply(CH.secondary[[t]],
                            1,
                            function(x){length(which(is.finite(x)))})
   }
   
   
   # NOTE - this only worked for Navajo and Grand Canyon
   # this makes a tibble of 0 x 6 when running on Zuni CH
   obs.dat <- obs.dat[-which(is.na(obs.dat$State)), ]
   
   
 # save Navajo data
   save.image("Data/Updated/Navajo_CaptureHistories_Jan2020.RData")
   
   
## Format data for Grand Canyon --------------------------------------------- ##
## -------------------------------------------------------------------------- ##
   
   
 # clear environment and reload necessary files
 # see above for comments
   
   
   CH.secondary <- multisite.CJS.capture.history.fun(
     dirty.data   = southwest.dirty, 
     cleaned.data = southwest.final.clean,
     site.webs    = c("grandcanyon.e", "grandcanyon.m", "grandcanyon.t"), 
     species      = "PM")
   
   
   session.list <- multisite.session.list.fun(
     dirty.data = southwest.dirty,
     site.webs  = c("grandcanyon.e", "grandcanyon.m", "grandcanyon.t")
   )
   

   CH.primary <- primary.ch.fun(CH.secondary)
   

   for (m in 1:length(CH.secondary)) {
     colnames(CH.secondary[[m]]) <- 1:dim(CH.secondary[[m]])[2]
   }
   

   n.sec.occ <- unlist(lapply(CH.secondary, function(x) {
     dim(x)[2]
   }))
   

   covariate.data <- multisite.monthly.covariate.fun(
     cleaned.data = southwest.final.clean, # cleaned data
     CH.secondary = CH.secondary,          # as monthly list
     tags = rownames(CH.secondary[[1]]),
     by.sitetag = TRUE,
     sessions = session.list$all.sessions,
     temporal.data = sw.temp.data,
     site.webs = c("grandcanyon.e", "grandcanyon.m", "grandcanyon.t"),
     cov.list = list(prcp = 0, prcp3 = 0, prcp6 = 0, prcp12 = 0, #prcp18 = 0,
                     ndvi = 0, ndvi3 = 0, ndvi6 = 0, ndvi12 = 0, #ndvi18 = 0,
                     temp = 0, temp3 = 0, temp6 = 0, temp12 = 0, #temp18 = 0,
                     tmin = 0, tmin3 = 0, tmin6 = 0, tmin12 = 0, #tmin18 = 0,
                     tmax = 0, tmax3 = 0, tmax6 = 0, tmax12 = 0, #tmax18 = 0,
                     swe  = 0, swe3  = 0, swe6  = 0, swe12  = 0,
                     swewinter = 0) #swe18  = 0)
   )
   

   p.or.c <- multisite.p.or.c.array.fun(CH.secondary, 
                                        covariate.data$temporal.covariates)
   
   
## Convert to long format --------------------------------------------------- ##
   
   
   obs.dat <- monthly.longdata.CH.fun(CH.secondary, 
                                      covariate.data$temporal.covariates, 
                                      covariate.data$individual.covariates)
   

   monthlyCH <- monthly.primary.CH.fun(CH.primary,
                                       covariate.data$temporal.covariates)
   

   webmonths <- list()
   

   for (i in 1:(length(session.list) - 1)) {
     x <- match(session.list[[i + 1]], 
                covariate.data$temporal.covariates$session)
     webmonths[[i]] <- x[which(is.finite(x))]
   }
   
   names(webmonths) <- names(session.list)[-1]
   

   months.trapped.mat <- matrix(NA, 
                                nrow = dim(CH.secondary[[1]])[1], 
                                ncol = max(unlist(lapply(webmonths, length))))
   

   length.months.trapped <- numeric()
   

   for (i in 1:dim(months.trapped.mat)[1]) {
     webnam <- covariate.data$individual.covariates$web[i] 
     webi   <- which(names(webmonths) == paste("web", webnam, sep = "."))
     length.months.trapped[i] <- length(webmonths[[webi]])
     months.trapped.mat[i, 1:length.months.trapped[i]] <- webmonths[[webi]]
   }
 
   
## Arrange covariates in an array ------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
 # create array for prcp
   prcp.covariate.array <- array(c(
     covariate.data$prcp_0, # 1
     covariate.data$prcp3_0, # 2
     covariate.data$prcp6_0, # 3
     covariate.data$prcp12_0 # 4
     #covariate.data$prcp18_0 # 5
   ), dim = c(dim(covariate.data$prcp_0), 4)) #5))
   
   
 # create array for ndvi
   ndvi.covariate.array <- array(c(
     covariate.data$ndvi_0, # 1
     covariate.data$ndvi3_0, # 2
     covariate.data$ndvi6_0, # 3
     covariate.data$ndvi12_0 # 4
     #covariate.data$ndvi18_0 # 5
   ), dim = c(dim(covariate.data$ndvi_0), 4)) #5))
   
   
 # create array for temp
   temp.covariate.array <- array(c(
     covariate.data$temp_0, # 1
     covariate.data$temp3_0, # 2
     covariate.data$temp6_0, # 3
     covariate.data$temp12_0 # 4
     #covariate.data$temp18_0 # 5
   ), dim = c(dim(covariate.data$temp_0), 4)) #5))
   
   
 # create array for tmin
   tmin.covariate.array <- array(c(
     covariate.data$tmin_0, # 1
     covariate.data$tmin3_0, # 2
     covariate.data$tmin6_0, # 3
     covariate.data$tmin12_0 # 4
     #covariate.data$tmin18_0 # 5
   ), dim = c(dim(covariate.data$tmin_0), 4)) #5))
   
   
 # create array for tmax
   tmax.covariate.array <- array(c(
     covariate.data$tmax_0, # 1
     covariate.data$tmax3_0, # 2
     covariate.data$tmax6_0, # 3
     covariate.data$tmax12_0 # 4
     #covariate.data$tmax18_0 # 5
   ), dim = c(dim(covariate.data$tmax_0), 4)) #5))
   
   
 # create array for swe
   swe.covariate.array <- array(c(
     covariate.data$swe_0, # 1
     covariate.data$swe3_0, # 2
     covariate.data$swe6_0, # 3
     covariate.data$swe12_0, # 4
     ##covariate.data$swe18_0 # 5
     covariate.data$swewinter_0
   ), dim = c(dim(covariate.data$swe_0), 5)) #5))
   
   
   # take care of problem of differing number of secondary occasions
   n.sec.occ <- matrix(NA,
                       nrow = dim(CH.primary)[1],
                       ncol = dim(CH.primary)[2])
   
   for(t in 1:dim(n.sec.occ)[2]) {
     n.sec.occ[,t] <- apply(CH.secondary[[t]],
                            1,
                            function(x){length(which(is.finite(x)))})
   }
   
   
   # NOTE - this only worked for Navajo and Grand Canyon
   # this makes a tibble of 0 x 6 when running on Zuni CH
   obs.dat <- obs.dat[-which(is.na(obs.dat$State)), ]
   
   
 # save Grand Canyon data
   save.image("Data/Updated/GrandCanyon_CaptureHistories_Jan2020.RData")  
   
   
## -------------------------------------------------------------------------- ##