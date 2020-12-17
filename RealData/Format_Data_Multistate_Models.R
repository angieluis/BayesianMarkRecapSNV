##################################################################
## Code to format data for multistate infection model
## for multiple sites in the same model
##################################################################




### right now these are the NM and AZ site/webs that we have environmental data for:
all.site.webs <- c("grandcanyon.e", "grandcanyon.m", "grandcanyon.t", 
                   "limestone.c1", "limestone.s1", "limestone.s2",  
                   "navajo.1",   "navajo.2",      
                   "placitas.1", "placitas.2",  "placitas.3",     
                   "walnutcreek.a", "walnutcreek.w", 
                   "zuni.1",  "zuni.2")    

setwd("~/Documents/JAGS/BayesianMarkRecapSNV/RealData")

load("AllCaptureData.RData")

source("01_sw_data_functions_more.R")

# read in normalized covariates
sw.temp.data <- read.csv("updated_southwest_covariates_norm.csv")

# load diversity/density data
# See "DiversityListCode.R"
load("DiversityLongdataInterpolated.RData")




# Format data and save for each site separately and then 
# will combine with combine data functions
  
# Grand Canyon -----------------------------------------------------------------------------##

ms.CH.secondary_GC <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("grandcanyon.e","grandcanyon.m","grandcanyon.t"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_GC <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs  = c("grandcanyon.e","grandcanyon.m","grandcanyon.t")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_GC <- primary.MSch.fun(ms.CH.secondary_GC)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_GC)) {
  colnames(ms.CH.secondary_GC[[m]]) <- 1:dim(ms.CH.secondary_GC[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_GC <- unlist(lapply(ms.CH.secondary_GC, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_GC <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_GC,          # as monthly list
    # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_GC[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_GC$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("grandcanyon.e","grandcanyon.m","grandcanyon.t"),
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
p.or.c_GC <- multisite.p.or.c.array.fun(ms.CH.secondary_GC, 
                                     covariate.data_GC$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_GC <- monthly.longdata.CH.fun(ms.CH.secondary_GC, 
                                   covariate.data_GC$temporal.covariates, 
                                   covariate.data_GC$individual.covariates)


monthlyCH_GC <- monthly.primary.CH.fun(ms.CH.primary_GC,
                                    covariate.data_GC$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_GC <- list()
# loop through sessions   
for (i in 1:(length(session.list_GC) - 1)) {
  x <- match(session.list_GC[[i + 1]], 
             covariate.data_GC$temporal.covariates$session)
  webmonths_GC[[i]] <- x[which(is.finite(x))]
}

names(webmonths_GC) <- names(session.list_GC)[-1]


# matrix for months trapped
months.trapped.mat_GC <- matrix(NA, 
                             nrow = dim(ms.CH.secondary_GC[[1]])[1], 
                             ncol = max(unlist(lapply(webmonths_GC, length))))


# create numeric
length.months.trapped_GC <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_GC)[1]) {
  # this is a factor currently
  webnam <- covariate.data_GC$individual.covariates$web[i] 
  webi   <- which(names(webmonths_GC) == paste("web", webnam, sep = "."))
  length.months.trapped_GC[i] <- length(webmonths_GC[[webi]])
  months.trapped.mat_GC[i, 1:length.months.trapped_GC[i]] <- webmonths_GC[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_GC <- matrix(NA,
                    nrow = dim(ms.CH.primary_GC)[1],
                    ncol = dim(ms.CH.primary_GC)[2])

for(t in 1:dim(n.sec.occ_GC)[2]) {
  n.sec.occ_GC[,t] <- apply(ms.CH.secondary_GC[[t]],
                         1,
                         function(x){length(which(is.finite(x)))})
}





# Limestone  -----------------------------------------------------------------------------##

ms.CH.secondary_L <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("limestone.c1", "limestone.s1", "limestone.s2"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_L <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs = c("limestone.c1", "limestone.s1", "limestone.s2")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_L <- primary.MSch.fun(ms.CH.secondary_L)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_L)) {
  colnames(ms.CH.secondary_L[[m]]) <- 1:dim(ms.CH.secondary_L[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_L <- unlist(lapply(ms.CH.secondary_L, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_L <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_L,          # as monthly list
  # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_L[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_L$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("limestone.c1", "limestone.s1", "limestone.s2"),
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
p.or.c_L <- multisite.p.or.c.array.fun(ms.CH.secondary_L, 
                                        covariate.data_L$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_L <- monthly.longdata.CH.fun(ms.CH.secondary_L, 
                                      covariate.data_L$temporal.covariates, 
                                      covariate.data_L$individual.covariates)


monthlyCH_L <- monthly.primary.CH.fun(ms.CH.primary_L,
                                       covariate.data_L$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_L <- list()
# loop through sessions   
for (i in 1:(length(session.list_L) - 1)) {
  x <- match(session.list_L[[i + 1]], 
             covariate.data_L$temporal.covariates$session)
  webmonths_L[[i]] <- x[which(is.finite(x))]
}

names(webmonths_L) <- names(session.list_L)[-1]


# matrix for months trapped
months.trapped.mat_L <- matrix(NA, 
                                nrow = dim(ms.CH.secondary_L[[1]])[1], 
                                ncol = max(unlist(lapply(webmonths_L, length))))


# create numeric
length.months.trapped_L <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_L)[1]) {
  # this is a factor currently
  webnam <- covariate.data_L$individual.covariates$web[i] 
  webi   <- which(names(webmonths_L) == paste("web", webnam, sep = "."))
  length.months.trapped_L[i] <- length(webmonths_L[[webi]])
  months.trapped.mat_L[i, 1:length.months.trapped_L[i]] <- webmonths_L[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_L <- matrix(NA,
                       nrow = dim(ms.CH.primary_L)[1],
                       ncol = dim(ms.CH.primary_L)[2])

for(t in 1:dim(n.sec.occ_L)[2]) {
  n.sec.occ_L[,t] <- apply(ms.CH.secondary_L[[t]],
                            1,
                            function(x){length(which(is.finite(x)))})
}




# Navajo  -----------------------------------------------------------------------------##

ms.CH.secondary_N <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("navajo.1",   "navajo.2"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_N <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs = c("navajo.1",   "navajo.2")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_N <- primary.MSch.fun(ms.CH.secondary_N)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_N)) {
  colnames(ms.CH.secondary_N[[m]]) <- 1:dim(ms.CH.secondary_N[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_N <- unlist(lapply(ms.CH.secondary_N, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_N <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_N,          # as monthly list
  # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_N[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_N$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("navajo.1",   "navajo.2"),
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
p.or.c_N <- multisite.p.or.c.array.fun(ms.CH.secondary_N, 
                                       covariate.data_N$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_N <- monthly.longdata.CH.fun(ms.CH.secondary_N, 
                                     covariate.data_N$temporal.covariates, 
                                     covariate.data_N$individual.covariates)


monthlyCH_N <- monthly.primary.CH.fun(ms.CH.primary_N,
                                      covariate.data_N$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_N <- list()
# loop through sessions   
for (i in 1:(length(session.list_N) - 1)) {
  x <- match(session.list_N[[i + 1]], 
             covariate.data_N$temporal.covariates$session)
  webmonths_N[[i]] <- x[which(is.finite(x))]
}

names(webmonths_N) <- names(session.list_N)[-1]


# matrix for months trapped
months.trapped.mat_N <- matrix(NA, 
                               nrow = dim(ms.CH.secondary_N[[1]])[1], 
                               ncol = max(unlist(lapply(webmonths_N, length))))


# create numeric
length.months.trapped_N <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_N)[1]) {
  # this is a factor currently
  webnam <- covariate.data_N$individual.covariates$web[i] 
  webi   <- which(names(webmonths_N) == paste("web", webnam, sep = "."))
  length.months.trapped_N[i] <- length(webmonths_N[[webi]])
  months.trapped.mat_N[i, 1:length.months.trapped_N[i]] <- webmonths_N[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_N <- matrix(NA,
                      nrow = dim(ms.CH.primary_N)[1],
                      ncol = dim(ms.CH.primary_N)[2])

for(t in 1:dim(n.sec.occ_N)[2]) {
  n.sec.occ_N[,t] <- apply(ms.CH.secondary_N[[t]],
                           1,
                           function(x){length(which(is.finite(x)))})
}




# Placitas  -----------------------------------------------------------------------------##

ms.CH.secondary_P <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("placitas.1", "placitas.2",  "placitas.3"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_P <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs = c("placitas.1", "placitas.2",  "placitas.3")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_P <- primary.MSch.fun(ms.CH.secondary_P)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_P)) {
  colnames(ms.CH.secondary_P[[m]]) <- 1:dim(ms.CH.secondary_P[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_P <- unlist(lapply(ms.CH.secondary_P, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_P <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_P,          # as monthly list
  # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_P[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_P$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("placitas.1", "placitas.2",  "placitas.3"),
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
p.or.c_P <- multisite.p.or.c.array.fun(ms.CH.secondary_P, 
                                       covariate.data_P$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_P <- monthly.longdata.CH.fun(ms.CH.secondary_P, 
                                     covariate.data_P$temporal.covariates, 
                                     covariate.data_P$individual.covariates)


monthlyCH_P <- monthly.primary.CH.fun(ms.CH.primary_P,
                                      covariate.data_P$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_P <- list()
# loop through sessions   
for (i in 1:(length(session.list_P) - 1)) {
  x <- match(session.list_P[[i + 1]], 
             covariate.data_P$temporal.covariates$session)
  webmonths_P[[i]] <- x[which(is.finite(x))]
}

names(webmonths_P) <- names(session.list_P)[-1]


# matrix for months trapped
months.trapped.mat_P <- matrix(NA, 
                               nrow = dim(ms.CH.secondary_P[[1]])[1], 
                               ncol = max(unlist(lapply(webmonths_P, length))))


# create numeric
length.months.trapped_P <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_P)[1]) {
  # this is a factor currently
  webnam <- covariate.data_P$individual.covariates$web[i] 
  webi   <- which(names(webmonths_P) == paste("web", webnam, sep = "."))
  length.months.trapped_P[i] <- length(webmonths_P[[webi]])
  months.trapped.mat_P[i, 1:length.months.trapped_P[i]] <- webmonths_P[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_P <- matrix(NA,
                      nrow = dim(ms.CH.primary_P)[1],
                      ncol = dim(ms.CH.primary_P)[2])

for(t in 1:dim(n.sec.occ_P)[2]) {
  n.sec.occ_P[,t] <- apply(ms.CH.secondary_P[[t]],
                           1,
                           function(x){length(which(is.finite(x)))})
}




# Walnut Creek  -----------------------------------------------------------------------------##

ms.CH.secondary_WC <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("walnutcreek.a", "walnutcreek.w"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_WC <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs = c("walnutcreek.a", "walnutcreek.w")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_WC <- primary.MSch.fun(ms.CH.secondary_WC)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_WC)) {
  colnames(ms.CH.secondary_WC[[m]]) <- 1:dim(ms.CH.secondary_WC[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_WC <- unlist(lapply(ms.CH.secondary_WC, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_WC <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_WC,          # as monthly list
  # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_WC[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_WC$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("walnutcreek.a", "walnutcreek.w"),
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
p.or.c_WC <- multisite.p.or.c.array.fun(ms.CH.secondary_WC, 
                                       covariate.data_WC$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_WC <- monthly.longdata.CH.fun(ms.CH.secondary_WC, 
                                     covariate.data_WC$temporal.covariates, 
                                     covariate.data_WC$individual.covariates)


monthlyCH_WC <- monthly.primary.CH.fun(ms.CH.primary_WC,
                                      covariate.data_WC$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_WC <- list()
# loop through sessions   
for (i in 1:(length(session.list_WC) - 1)) {
  x <- match(session.list_WC[[i + 1]], 
             covariate.data_WC$temporal.covariates$session)
  webmonths_WC[[i]] <- x[which(is.finite(x))]
}

names(webmonths_WC) <- names(session.list_WC)[-1]


# matrix for months trapped
months.trapped.mat_WC <- matrix(NA, 
                               nrow = dim(ms.CH.secondary_WC[[1]])[1], 
                               ncol = max(unlist(lapply(webmonths_WC, length))))


# create numeric
length.months.trapped_WC <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_WC)[1]) {
  # this is a factor currently
  webnam <- covariate.data_WC$individual.covariates$web[i] 
  webi   <- which(names(webmonths_WC) == paste("web", webnam, sep = "."))
  length.months.trapped_WC[i] <- length(webmonths_WC[[webi]])
  months.trapped.mat_WC[i, 1:length.months.trapped_WC[i]] <- webmonths_WC[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_WC <- matrix(NA,
                      nrow = dim(ms.CH.primary_WC)[1],
                      ncol = dim(ms.CH.primary_WC)[2])

for(t in 1:dim(n.sec.occ_WC)[2]) {
  n.sec.occ_WC[,t] <- apply(ms.CH.secondary_WC[[t]],
                           1,
                           function(x){length(which(is.finite(x)))})
}




# Zuni  -----------------------------------------------------------------------------##

ms.CH.secondary_Z <- multisite.MS.capture.history.fun(
  dirty.data = southwest.dirty,
  cleaned.data = southwest.final.clean, # 'cleaned'
  site.webs = c("zuni.1",  "zuni.2"),
  species="PM",
  SNV.unknown.state=FALSE # if want -9 for SNV it's own state ("3"), otherwise assume unknowns are SNV negative ("1")
)

# make a list of every session trapped from dirty data
# list of 3 vectors - length of web plus one
# the first element of the list is all the sessions for all of the sites
# all other elements are sessions trapped for that specific web
session.list_Z <- multisite.session.list.fun(
  dirty.data = southwest.dirty,
  site.webs = c("zuni.1",  "zuni.2")
)

# makes a primary capture history from secondary occasions
ms.CH.primary_Z <- primary.MSch.fun(ms.CH.secondary_Z)

# rename secondary occasions from 1 to number of days
# looking at column names and changing them
for (m in 1:length(ms.CH.secondary_Z)) {
  colnames(ms.CH.secondary_Z[[m]]) <- 1:dim(ms.CH.secondary_Z[[m]])[2]
} 


# number of secondary occasions - maximum of secondary occasions across webs
n.sec.occ_Z <- unlist(lapply(ms.CH.secondary_Z, function(x) {
  dim(x)[2]
}))



# apply covariate function to normalized data
covariate.data_Z <- multisite.monthly.covariate.fun(
  cleaned.data = southwest.final.clean, # cleaned data
  CH.secondary = ms.CH.secondary_Z,          # as monthly list
  # in CH.secondary row names are tags, tag names line up to CH.secondary
  tags = rownames(ms.CH.secondary_Z[[1]]),
  by.sitetag = TRUE, 
  sessions = session.list_Z$all.sessions,
  temporal.data = sw.temp.data, 
  multistate=TRUE,
  diversity.data = scaled.MNAs.diversity.interp.longdata,
  remove.na=TRUE,
  site.webs = c("zuni.1",  "zuni.2"),
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
p.or.c_Z <- multisite.p.or.c.array.fun(ms.CH.secondary_Z, 
                                        covariate.data_Z$temporal.covariates)

## Convert to long format --------------------------------------------------- ##

obs.dat_Z <- monthly.longdata.CH.fun(ms.CH.secondary_Z, 
                                      covariate.data_Z$temporal.covariates, 
                                      covariate.data_Z$individual.covariates)


monthlyCH_Z <- monthly.primary.CH.fun(ms.CH.primary_Z,
                                       covariate.data_Z$temporal.covariates)


# accounts for how certain webs weren't trapped certain months, that are
# within the long data set (i.e., shit wasn't trapped at the same time - month)
# we use all the covariates but we only trapped some of those times, need
# covariates for the observations process, simulating state process over all
# months, even those not trapped
webmonths_Z <- list()
# loop through sessions   
for (i in 1:(length(session.list_Z) - 1)) {
  x <- match(session.list_Z[[i + 1]], 
             covariate.data_Z$temporal.covariates$session)
  webmonths_Z[[i]] <- x[which(is.finite(x))]
}

names(webmonths_Z) <- names(session.list_Z)[-1]


# matrix for months trapped
months.trapped.mat_Z <- matrix(NA, 
                                nrow = dim(ms.CH.secondary_Z[[1]])[1], 
                                ncol = max(unlist(lapply(webmonths_Z, length))))


# create numeric
length.months.trapped_Z <- numeric()
# loop through months trapped
for (i in 1:dim(months.trapped.mat_Z)[1]) {
  # this is a factor currently
  webnam <- covariate.data_Z$individual.covariates$web[i] 
  webi   <- which(names(webmonths_Z) == paste("web", webnam, sep = "."))
  length.months.trapped_Z[i] <- length(webmonths_Z[[webi]])
  months.trapped.mat_Z[i, 1:length.months.trapped_Z[i]] <- webmonths_Z[[webi]]
}

# take care of problem of differing number of secondary occasions
n.sec.occ_Z <- matrix(NA,
                       nrow = dim(ms.CH.primary_Z)[1],
                       ncol = dim(ms.CH.primary_Z)[2])

for(t in 1:dim(n.sec.occ_Z)[2]) {
  n.sec.occ_Z[,t] <- apply(ms.CH.secondary_Z[[t]],
                            1,
                            function(x){length(which(is.finite(x)))})
}






#-------------------------------------------------------------------###
 # Now combine all the sites together

obs.dat <- combine.obsdat.fun(obs.dats = list(obs.dat_GC,obs.dat_L,obs.dat_N,obs.dat_P,obs.dat_WC,obs.dat_Z), individual.covariates = list(covariate.data_GC$individual.covariates,covariate.data_L$individual.covariates,covariate.data_N$individual.covariates, covariate.data_P$individual.covariates,covariate.data_WC$individual.covariates,covariate.data_Z$individual.covariates))

## remove NA rows of obs.dat (when not trapped that day)
obs.dat <- obs.dat %>%
  dplyr::filter(!is.na(State))

## replace State 0 (not observed) to 3 for multistate infection model
obs.dat$State <- replace(obs.dat$State,obs.dat$State==0,3)

covariate.data <- combine.covariates.fun(covariate.data = list(covariate.data_GC,covariate.data_L,covariate.data_N,covariate.data_P,covariate.data_WC, covariate.data_Z))

monthlyCH <- combine.monthlyCH.fun(monthlyCH = list(monthlyCH_GC,monthlyCH_L,monthlyCH_N,monthlyCH_P,monthlyCH_WC,monthlyCH_Z))

p.or.c <- combine.porc.fun(p.or.c = list(p.or.c_GC, p.or.c_L, p.or.c_N,p.or.c_P,p.or.c_WC, p.or.c_Z))

n.sec.occ <- combine.n.sec.occ.func(n.sec.occ = list(n.sec.occ_GC,n.sec.occ_L,n.sec.occ_N,n.sec.occ_P,n.sec.occ_WC,n.sec.occ_Z))

months.trapped.mat <- combine.months.trapped.func(months.trapped.mat = list(months.trapped.mat_GC,months.trapped.mat_L, months.trapped.mat_N,months.trapped.mat_P,months.trapped.mat_WC,months.trapped.mat_Z))

length.months.trapped <- c(length.months.trapped_GC,length.months.trapped_L,length.months.trapped_N,length.months.trapped_P,length.months.trapped_WC,length.months.trapped_Z)


save(obs.dat, covariate.data, monthlyCH, p.or.c, n.sec.occ, months.trapped.mat, length.months.trapped,file="CombinedAllSitesMSData.RData")

  