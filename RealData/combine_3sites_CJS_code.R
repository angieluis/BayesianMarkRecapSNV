#################################################################################################
#  Use functions in "combine_data_functions.R" to combine data for Grand Canyon, Navajo, and Zuni
# This is all for CJS models
#################################################################################################


load("ThreeSiteCJSData.RData")
### Data generated from previous functions: covariates.data list, p.or.c arrays, etc for the 3 sites

source("Code/combine_data_functions.R") # or now in "01_sw_data_functions_more.R"


obs.dat <- combine.obsdat.fun(obs.dats = list(obs.dat_GC,obs.dat_N,obs.dat_Z), individual.covariates = list(covariate.data_GC$individual.covariates,covariate.data_N$individual.covariates, covariate.data_Z$individual.covariates))

covariate.data <- combine.covariates.fun(covariate.data = list(covariate.data_GC,covariate.data_N, covariate.data_Z))

monthlyCH <- combine.monthlyCH.fun(monthlyCH = list(monthlyCH_GC,monthlyCH_N,monthlyCH_Z))

p.or.c <- combine.porc.fun(p.or.c = list(p.or.c_GC, p.or.c_N, p.or.c_Z))

n.sec.occ <- combine.n.sec.occ.func(n.sec.occ = list(n.sec.occ_GC,n.sec.occ_N,n.sec.occ_Z))

months.trapped.mat <- combine.months.trapped.func(months.trapped.mat = list(months.trapped.mat_GC,months.trapped.mat_N,months.trapped.mat_Z))

length.months.trapped <- c(length.months.trapped_GC,length.months.trapped_N,length.months.trapped_Z)


save(obs.dat, covariate.data, monthlyCH, p.or.c, n.sec.occ, months.trapped.mat, length.months.trapped,file="CombinedThreeSiteCJSData.RData")



  

  
