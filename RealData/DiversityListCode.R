
# Create a list of MNA, MNI, and diversity indices for all sites/webs------------------------##

load("AllCaptureData.RData")

source("01_sw_data_functions_more.R")



MNAs.diversity.list <- list()
site.webs <- sort(unique(southwest.final.clean$site_web))
site.webs <- sub("_",".",site.webs)
for(i in 1:length(site.webs)){
  site <- unlist(strsplit(site.webs[i],"[.]"))[1]
  web <- unlist(strsplit(site.webs[i],"[.]"))[2]
  MNAs.diversity.list[[i]] <- diversity.df.function(
    data = southwest.final.clean,   
    site = site,  
    web = web,  
    sessions = multisite.session.list.fun(dirty.data=southwest.dirty,
                                          site.webs=site.webs[i])[[1]],  
    interpolate=FALSE, 
    scale=FALSE, 
    include.pm = TRUE)
}
names(MNAs.diversity.list) <- site.webs
save(MNAs.diversity.list,file="DiversityDataList.RData")

#### Make it long data to match sw.temp.data

MNAs.diversity.longdata <- MNAs.diversity.list[[1]]
for(i in 2:length(MNAs.diversity.list)){
  MNAs.diversity.longdata <- rbind(MNAs.diversity.longdata,MNAs.diversity.list[[i]])
}

## make max 1 for all covariates across all sites
scaled.MNAs.diversity.longdata <- MNAs.diversity.longdata[,1:7]  # assumes these columns refer to the times and the rest of the columns are the MNAs and diversities
for(i in 8:dim(MNAs.diversity.longdata)[2]){
  scaled.MNAs.diversity.longdata[,i] <- MNAs.diversity.longdata[,i]/max(MNAs.diversity.longdata[,i],na.rm=TRUE)
  names(scaled.MNAs.diversity.longdata)[i] <- names(MNAs.diversity.longdata)[i]
}

save(MNAs.diversity.longdata, scaled.MNAs.diversity.longdata, file="DiversityLongdata.RData")

###########################################################################################
## Same, but now interpolate. (Models can't have NAs)
###########################################################################################

# Create a list of MNA, MNI, and diversity indices for all sites/webs------------------------##

MNAs.diversity.list.interp <- list()
site.webs <- sort(unique(southwest.final.clean$site_web))
site.webs <- sub("_",".",site.webs)
for(i in 1:length(site.webs)){
  site <- unlist(strsplit(site.webs[i],"[.]"))[1]
  web <- unlist(strsplit(site.webs[i],"[.]"))[2]
  MNAs.diversity.list.interp[[i]] <- diversity.df.function(
    data = southwest.final.clean,   
    site = site,  
    web = web,  
    sessions = multisite.session.list.fun(dirty.data=southwest.dirty,
                                          site.webs=site.webs[i])[[1]],  
    interpolate=TRUE, 
    scale=FALSE, 
    include.pm = TRUE)
}
names(MNAs.diversity.list.interp) <- site.webs
save(MNAs.diversity.list.interp,file="DiversityDataListInterpolated.RData")

#### Make it long data to match sw.temp.data

MNAs.diversity.interp.longdata <- MNAs.diversity.list.interp[[1]]
for(i in 2:length(MNAs.diversity.list.interp)){
  MNAs.diversity.interp.longdata <- rbind(MNAs.diversity.interp.longdata,MNAs.diversity.list.interp[[i]])
}

## make max 1 for all covariates across all sites
scaled.MNAs.diversity.interp.longdata <- MNAs.diversity.interp.longdata[,1:7]  # assumes these columns refer to the times and the rest of the columns are the MNAs and diversities
for(i in 8:dim(MNAs.diversity.interp.longdata)[2]){
  scaled.MNAs.diversity.interp.longdata[,i] <- MNAs.diversity.interp.longdata[,i]/max(MNAs.diversity.interp.longdata[,i],na.rm=TRUE)
  names(scaled.MNAs.diversity.interp.longdata)[i] <- names(MNAs.diversity.interp.longdata)[i]
}

save(MNAs.diversity.interp.longdata, scaled.MNAs.diversity.interp.longdata, file="DiversityLongdataInterpolated.RData")
