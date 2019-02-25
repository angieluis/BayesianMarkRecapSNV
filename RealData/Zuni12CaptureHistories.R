### Get dates from the uncleaned data (NAs not removed)
head(UNMcaptures)
UNMcaptures$date <- as.Date(gsub(" ", "", UNMcaptures$Date),format="%m/%d/%Y")
sort(unique(UNMcaptures$date))
Zuni.dates <- sort(unique(UNMcaptures$date[which(UNMcaptures$site=="Zuni")]))
Zuni.sessions <- sort(unique(UNMcaptures$Session[which(UNMcaptures$site=="Zuni")]))

time.int <- diff(Zuni.dates)
first.dates <- Zuni.dates[c(1,1+which(time.int>1))]
length(first.dates)

primary.time.int.weeks <- diff(first.dates)/7
Zuni.primary.time.int.weeks <- primary.time.int.weeks

### use those dates on the cleaned data:
UNMdata$date <- as.Date(gsub(" ", "", UNMdata$Date),format="%m/%d/%Y")
Zuni.data <- UNMdata[which(UNMdata$site=="Zuni"),]
Zuni12.data <- Zuni.data[-which(Zuni.data$web==3),]

# I want to separate out deermice. Since some of the repeaters might have -9 under letter_2, I will go by tag. These are the tag numbers I want to keep (are pema).
Zuni12.pema.tags <- sort(unique(Zuni12.data$tag[which(Zuni12.data$letter_2=="PM")]))
# a couple have two tags and one is seen elsewhere. make same id?
Zuni12.pema.tags <- Zuni12.pema.tags[-which(Zuni12.pema.tags==-9)]

pema.ind <- numeric()
for(i in 1:length(Zuni12.pema.tags)){
  pema.ind <- c(pema.ind,which(Zuni12.data$tag==Zuni12.pema.tags[i]))
}
Zuni12.pema.data <- Zuni12.data[pema.ind,]



Zuni12.Session.days <- list()
for(i in 1:length(first.dates)){
  Zuni12.Session.days[[i]] <- Zuni.dates[which(Zuni.dates==first.dates[i]):ifelse(i==length(first.dates),length(Zuni.dates),(which(Zuni.dates==first.dates[i+1])-1))]
}




################################################################################
########## Basic Robust Design CJS capture histories (0 or 1) AS LIST
# A matrix for each month (primary occasion), where column is day (secondary occasion) and row is individual. Each matrix has the same number of rows - all individuals have a row each month even though for most they will be zeros.
################################################################################
IDs <- sort (unique(Zuni12.pema.data$tag))
Session.days <- Zuni12.Session.days

Ch.list <- list()

for(m in 1:length(Session.days)){
  days <- Session.days[[m]]
  ch.mat <- matrix(NA,ncol=length(days),nrow=length(IDs))
  
  for(d in 1:length(days)){
    for(i in 1:length(IDs)){
      ch.mat[i,d] <- ifelse(length(which(Zuni12.pema.data$tag==IDs[i] & Zuni12.pema.data$date==days[d]))>0,1,0)
    }
  }
  dimnames(ch.mat) <- list(IDs,NULL)
  Ch.list[[m]] <- ch.mat
  cat("session = ", m, "\n")
}

Zuni12.pema.Ch.secondary <- Ch.list
  
  
Zuni12.primary.time.int.weeks <- Zuni.primary.time.int.weeks 
Zuni12.secondary.occasions <- unlist(lapply(Zuni12.Session.days,length))



####### Temporal Covariates 

#source("RobustCJSfunctions.R")
#Zuni12.monthly.temporal.covariates <- monthly.temporaldata.fun(capture.data=UNMcaptures, site="Zuni",web=c(1,2))

#Zuni12.weekly.temporal.covariates <- weekly.temporaldata.fun(dates=Zuni12.dates)
####### Individual Covariates

#Zuni12.individual.covariates <- individual.covariate.fun(Zuni12.pema.data, IDs, Zuni12.pema.Ch.secondary)


##############

save(Zuni12.pema.Ch.secondary,Zuni12.Session.days,Zuni12.primary.time.int.weeks, Zuni12.secondary.occasions, Zuni12.pema.data, Zuni.sessions, file="Zuni12pemaCH.RData")
