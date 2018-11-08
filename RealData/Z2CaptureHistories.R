### Get dates from the uncleaned data (NAs not removed)
head(UNMcaptures)
UNMcaptures$date <- as.Date(gsub(" ", "", UNMcaptures$Date),format="%m/%d/%Y")
sort(unique(UNMcaptures$date))
Z2.dates <- sort(unique(UNMcaptures$date[which(UNMcaptures$site=="Zuni"&UNMcaptures$web==2)]))


time.int <- diff(Z2.dates)
first.dates <- Z2.dates[c(1,1+which(time.int>1))]
length(first.dates)

primary.time.int.weeks <- diff(first.dates)/7
Z2.primary.time.int.weeks <- primary.time.int.weeks

### use those dates on the cleaned data:
UNMdata$date <- as.Date(gsub(" ", "", UNMdata$Date),format="%m/%d/%Y")
Z2.data <- UNMdata[which(UNMdata$site=="Zuni"&UNMdata$web==2),]

# I want to separate out deermice. Since some of the repeaters might have -9 under letter_2, I will go by tag. These are the tag numbers I want to keep (are pema).
Z2.pema.tags <- sort(unique(Z2.data$tag[which(Z2.data$letter_2=="PM")]))
# a couple have two tags and one is seen elsewhere. make same id?
Z2.pema.tags <- Z2.pema.tags[-which(Z2.pema.tags==-9)]

pema.ind <- numeric()
for(i in 1:length(Z2.pema.tags)){
  pema.ind <- c(pema.ind,which(Z2.data$tag==Z2.pema.tags[i]))
}
Z2.pema.data <- Z2.data[pema.ind,]



Session.days <- list()
for(i in 1:length(first.dates)){
  Session.days[[i]] <- Z2.dates[which(Z2.dates==first.dates[i]):ifelse(i==length(first.dates),length(Z2.dates),(which(Z2.dates==first.dates[i+1])-1))]
}




################################################################################
########## Basic Robust Design CJS capture histories (0 or 1)
# A matrix for each month (primary occasion), where column is day (secondary occasion) and row is individual. Each matrix has the same number of rows - all individuals have a row each month even though for most they will be zeros.
################################################################################
IDs <- sort (unique(Z2.pema.data$tag))

Ch.list <- list()

for(m in 1:length(Session.days)){
  days <- Session.days[[m]]
  ch.mat <- matrix(NA,ncol=length(days),nrow=length(IDs))
  
  for(d in 1:length(days)){
    for(i in 1:length(IDs)){
      ch.mat[i,d] <- ifelse(length(which(Z2.pema.data$tag==IDs[i] & Z2.pema.data$date==days[d]))>0,1,0)
    }
  }
  Ch.list[[m]] <- ch.mat
  cat("session = ", m, "\n")
}

Z2.pema.Ch.secondary <- Ch.list
  
  
####### Temporal Covariates 
# month
library(lubridate)

month <- factor(month(first.dates))

temporal.covariates <- data.frame(session=1:length(Session.days),month)
