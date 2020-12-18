
#dummy temporal covariate data frame:
temp.cov <- data.frame(covariate.prim=c(1,2,2,3), # 1 to number of sessions, some numbers are repeated if months not trapped (the second 2 wasn't trapped)
      month.session=1:4, #numbered from 1 to number of months time span between first occasion and last occasion
      covariate.month=1:4, # this is month of the year, repeats 1:12
      NDVI=c(.1,.2,.25,.3))
attach(temp.cov)
#categorical variable (factor)

########################################## month as a factor

# 1 prior for each month

for(u in 1:12){
  mean.phi[u] ~ dnorm(0, 0.4)T(-10,10)     # priors for monthly survival
}

months <- month.session[which(covariate.prim==m)] #if there are no gaps, months should be length 1, but if some months weren't trapped there may be multiple months' data which should be incorporated

denom <- 1 + exp(-mean.phi[covariate.month[which(covariate.prim==m)]])

phi[i,m] <- (1/prod(denom))^(1/length(months))



######### for continuous covariate (NDVI):
# define priors for a (intercept) and b (slope)

months <- month.session[which(covariate.prim==m)]

denom <- 1 + exp(-a - b*NDVI[months])

phi[i,m] <- (1/prod(denom))^(1/length(months))


##### for month and NDVI:
### I think the monthly factors could replace the intercept term (so get rid of a) (12 different intercepts (each month) + 1 slope for how NDVI affects survival)

months <- month.session[which(covariate.prim==m)]

denom <- 1 + exp(-mean.phi[covariate.month[which(covariate.prim==m)]] - b*NDVI[months])

phi[i,m] <- (1/prod(denom))^(1/length(months))



