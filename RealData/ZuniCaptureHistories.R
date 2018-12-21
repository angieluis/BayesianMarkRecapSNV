### Get dates from the uncleaned data (NAs not removed)
head(UNMcaptures)
UNMcaptures$date <- as.Date(gsub(" ", "", UNMcaptures$Date),format="%m/%d/%Y")
sort(unique(UNMcaptures$date))
Zuni.dates <- sort(unique(UNMcaptures$date[which(UNMcaptures$site=="Zuni")]))


time.int <- diff(Zuni.dates)
first.dates <- Zuni.dates[c(1,1+which(time.int>1))]
length(first.dates)

primary.time.int.weeks <- diff(first.dates)/7
Zuni.primary.time.int.weeks <- primary.time.int.weeks

### use those dates on the cleaned data:
UNMdata$date <- as.Date(gsub(" ", "", UNMdata$Date),format="%m/%d/%Y")
Zuni.data <- UNMdata[which(UNMdata$site=="Zuni"),]

# I want to separate out deermice. Since some of the repeaters might have -9 under letter_2, I will go by tag. These are the tag numbers I want to keep (are pema).
Zuni.pema.tags <- sort(unique(Zuni.data$tag[which(Zuni.data$letter_2=="PM")]))
# a couple have two tags and one is seen elsewhere. make same id?
Zuni.pema.tags <- Zuni.pema.tags[-which(Zuni.pema.tags==-9)]

pema.ind <- numeric()
for(i in 1:length(Zuni.pema.tags)){
  pema.ind <- c(pema.ind,which(Zuni.data$tag==Zuni.pema.tags[i]))
}
Zuni.pema.data <- Zuni.data[pema.ind,]



Zuni.Session.days <- list()
for(i in 1:length(first.dates)){
  Zuni.Session.days[[i]] <- Zuni.dates[which(Zuni.dates==first.dates[i]):ifelse(i==length(first.dates),length(Zuni.dates),(which(Zuni.dates==first.dates[i+1])-1))]
}




################################################################################
########## Basic Robust Design CJS capture histories (0 or 1)
# A matrix for each month (primary occasion), where column is day (secondary occasion) and row is individual. Each matrix has the same number of rows - all individuals have a row each month even though for most they will be zeros.
################################################################################
IDs <- sort (unique(Zuni.pema.data$tag))

Ch.list <- list()

for(m in 1:length(Session.days)){
  days <- Session.days[[m]]
  ch.mat <- matrix(NA,ncol=length(days),nrow=length(IDs))
  
  for(d in 1:length(days)){
    for(i in 1:length(IDs)){
      ch.mat[i,d] <- ifelse(length(which(Zuni.pema.data$tag==IDs[i] & Zuni.pema.data$date==days[d]))>0,1,0)
    }
  }
  dimnames(ch.mat) <- list(IDs,NULL)
  Ch.list[[m]] <- ch.mat
  cat("session = ", m, "\n")
}

Zuni.pema.Ch.secondary <- Ch.list
  
  
####### Temporal Covariates 

source("RobustCJSfunctions.R")
Zuni.monthly.temporal.covariates <- monthly.temporaldata.fun(data=UNMcaptures, site="Zuni")

Zuni.weekly.temporal.covariates <- weekly.temporaldata.fun(dates=Zuni.dates)
####### Individual Covariates

Zuni.individual.covariates <- individual.covariate.fun(Zuni.pema.data, IDs, Zuni.pema.Ch.secondary)


##############

save(Zuni.pema.Ch.secondary,Zuni.Session.days,Zuni.weekly.temporal.covariates,Zuni.monthly.temporal.covariates,Zuni.primary.time.int.weeks,Zuni.individual.covariates, file="ZunipemaCH.RData")
