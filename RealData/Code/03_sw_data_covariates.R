## Get and arrange covariates ----------------------------------------------- ##


 # source files
   source("Code/01_sw_data_functions.R")
   source("Code/02_sw_data_clean.R")
   
   
## Southwest NDVI ----------------------------------------------------------- ##


 # read in data (made earlier)
   sw <- read.csv("Data/Covariates/southwest_ndvi.csv", stringsAsFactors = FALSE)


 # get rid of superfluous information
   ndvi.sw <- sw %>%
     
     # just southwest sites 
     filter(!agent == "mtech") %>%
  
     # get rid of arizona (university of arizona) and sevilleta
     # no data/low densities of deer mice
     filter(!siteAB %in% c("sevilleta1", "sevilleta2", "sevilleta3","arizona")) %>%

     # everything to lowercase
     mutate_all(funs(tolower)) %>%
  
     # order columns
     dplyr::select(year, month, site, web, ndvi) %>%
  
     # arrange alphabetically
     arrange(site, web)


## Daymet for the southwest ------------------------------------------------- ##


 # get location data for all SW sites
   swsites <- download_daymet_batch(
      file_location = "Data/southwest_locations.csv",
      start         = 1993,
      end           = 2007,
      internal      = FALSE,
                                 
      # will store files in this folder                                
      path = "Data/DaymetBatchSW")


 # combine into one data frame with column for file name (aka site)
   swsites <- list.files(path       = "Data/Covariates/DaymetBatchSW",
                         pattern    = "*.csv",
                         full.names = T) %>%
     map_df(function(x) read_csv(x, 
                                 skip = 7, 
                                 col_types = cols(.default = "c")) %>%
              mutate(filename = gsub(".csv","",basename(x))))


 # make dataframe
   swday <- as.data.frame(swsites)


 # eliminate decimal from year and yday columns
   swday[, 1:2] <- sapply(swday[, 1:2], as.numeric)


 # convert day (yday) to month day year - note that yday is in decimal
   swday <- swday %>%
     mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))


 # eliminate certain columns (srad, vp, year, yday, dayl)
 # 7/27/2019 - updated to add swe as a covariate
 # 7/28/2019 - code updated to eliminate tmin and tmax, use mean temp instead
   swdaymet <- swday %>%
     
     # rename columns
     dplyr::rename(swe  = "swe (kg/m^2)", 
                   tmin = "tmin (deg c)",
                   tmax = "tmax (deg c)",
                   prcp = "prcp (mm/day)",
                   site = "filename") %>%
  
     # select certain columns
     # prcp and swe (mm), tmin and tmax (C)
     dplyr::select(date, site, prcp, swe, tmax, tmin)


 # change format of date...again (this will make group_by and summarize better)
   swdaymet$date       <- as.Date(swdaymet$date, format = "%Y-%m-%d")


 # split site into agent, site, and web
   swsep <- swdaymet %>%
  
     # separate the site column 
     separate(site, into = paste(c("agent", "site", "web", "start", "end"), 
                                 sep = "_")) %>%
  
     # order columns
     dplyr::select(date, site, web, prcp, swe, tmax, tmin) %>%
  
     # arrange alphabetically
     arrange(site, web) %>%
  
     # change to lower case
     mutate_all(tolower)


 # make prcp, swe, tmin, tmax numeric
   swsep[, 4:7] <- sapply(swsep[, 4:7], as.numeric)


 # make year month column
   swsep <- swsep %>%
     mutate(date = as.yearmon(date))


## Southwest precip --------------------------------------------------------- ##


 # sum monthly precipitation by site 
   swprcp <- swsep %>%
  
     # group by site, then web, then date
     group_by(site, web, date) %>%

     # sum precipitation by month
     summarise(sum_prcp = sum(prcp))
  
  
 # order columns so everything aligns in join
   swprcp <- swprcp %>%
  
     dplyr::select(date, site, web, sum_prcp) %>%
  
     # arrange alphabetically
     arrange(site, web)


## Southwest swe ------------------------------------------------------------ ##
   
   
 # max monthly swe by site
   swswe <- swsep %>%
      
      # group by site, then web, then date
      group_by(site, web, date) %>%
      
      # sum precipitation by month
      summarise(max_swe = max(swe))
   
   
 # order columns so everything aligns in join
   swswe <- swswe %>%
      
      dplyr::select(date, site, web, max_swe) %>%
      
      # arrange alphabetically
      arrange(site, web)
   
   
## Southwest tmin ----------------------------------------------------------- ##


 # average tmin by month
   swtmin <- swsep %>%
  
     # group by site, then web, then date
     group_by(site, web, date) %>%
  
     # sum precipitation by month
     summarise(mean_tmin = mean(tmin))


 # order columns so everything aligns in join
   swtmin <- swtmin %>%
  
     dplyr::select(date, site, web, mean_tmin) %>%
  
     # arrange alphabetically
     arrange(site, web)


## Southwest tmax ----------------------------------------------------------- ##


 # average monthly tmax by site and web
   swtmax <- swsep %>%
  
     # group by site, then web, then date
     group_by(site, web, date) %>%
  
     # sum precipitation by month
     summarise(mean_tmax = mean(tmax))


 # order columns so everything aligns in join
   swtmax <- swtmax %>%
  
     dplyr::select(date, site, web, mean_tmax) %>%
  
     # arrange alphabetically
     arrange(site, web)


## Combine dataframes ------------------------------------------------------- ##


 # combine year and month
   ndvi.sw <- ndvi.sw %>%
      
      # combine
      mutate(date = as.yearmon(paste(year, " ", month), format = "%Y %m"))


 # rename, drop, and reorder
   swndvi <- ndvi.sw %>%
      
      # drop and reorder
      dplyr::select(date, site, web, ndvi)


 # merge the dataframes based on yearmon, site, and web
   southwest.merge <- swtmin %>%
      left_join(swtmax) %>%
      left_join(swprcp) %>%
      left_join(swndvi) %>%
      left_join(swswe)
   

## Interpolate NDVI NAs ----------------------------------------------------- ##


 # change ndvi from character to numeric but ignore NAs
   southwest.merge$ndvi <- as.numeric(southwest.merge$ndvi)
   
 
 # compute ndvi NAs based on monthly average throughout the years drop 2007
   sw.merge <- southwest.merge %>%
     filter(!grepl("2007", date)) %>%
     group_by(date) %>%
     mutate(new_ndvi = replace(ndvi, is.na(ndvi), mean(ndvi, na.rm = TRUE))) %>%
     dplyr::select(-ndvi)
   
   
 # rename columns
   sw.merge <- sw.merge %>% 
      dplyr::rename(tmin = mean_tmin,
                    tmax = mean_tmax,
                    prcp = sum_prcp,
                    swe  = max_swe,
                    ndvi = new_ndvi) 
   
   
## Calculate sum of SWE for spring months ----------------------------------- ##

   
 # does it work for everything???
   sw.swewinter <- sw.merge %>%
      
      # mutate a new column for month 1-12
      mutate(month = month(date)) %>%
      
      # group by site and web so it calculates appropriate value for each web
      group_by(site, web) %>%
      
      # new column
      mutate(swewinter = 
                
                # sum snow from oct to mar
                if_else(month == 3, 
                        rollsum(swe, 6, align = "right", fill = NA),
                        
                        # sum snow from oct to apr
                        if_else(month == 4, 
                                rollsum(swe, 7, align = "right", fill = NA),
                                
                                # sum snow from oct to may
                                if_else(month == 5, 
                                        rollsum(swe, 8, align = "right", fill = NA),
                                        
                                        # sum snow from oct to jun       
                                        if_else(month == 6, 
                                                rollsum(swe, 9, align = "right", fill = NA),
                                                
                                                # sum snow from oct to jul
                                                if_else(month == 7, 
                                                        rollsum(swe, 10, align = "right", fill = NA),
                                                        
                                                        # sum snow from oct to aug
                                                        if_else(month == 8, 
                                                                rollsum(swe, 11, align = "right", fill = NA),
                                                                
                                                                # sum snow from oct to sep         
                                                                if_else(month == 9, 
                                                                        rollsum(swe, 12, align = "right", fill = NA),
                                                                        
                                                                        # put 0 for everything else
                                                                        0)))))))) %>%
      
      # organize dataframe
      select(date, month, site, web, tmin, tmax, prcp, ndvi, swe, swewinter)
                                
   
## Calculate mean temp from min and max ------------------------------------- ##
   
 
 # ask Solomon if there are issues with this methodology
 # 8/29/2019 - he says nope but max and mins are probably more meaningful
 # per communication
   sw.temp.mean <- sw.swewinter %>%
      mutate(temp = rowMeans(cbind(tmin, tmax)))
   
   
## Seasons ------------------------------------------------------------------ ##
   
   
 # seasons as categories (i.e., 1 through 4 corresponding to winter, spring, 
 # summer, and fall, respectively)
 # angie coded seasons into the monthly covariate function
 # 1 = december - march: accounts for ENSO events
 # 2 = april - june
 # 3 = july - september: accounts for monsoon events
 # 4 = october - november
 # angie tackled seasons in her code

   
## Save to file before normalizing ------------------------------------------ ##
 

 # rename
   southwest <- sw.temp.mean
   

 # write to csv for later use in covariate function - NOT NORMALIZED
   write.csv(southwest, "Data/southwest_covariates.csv")

  
## Time lags ---------------------------------------------------------------- ##
    

 # what i have
      # average ndvi (values 0-1, 180 m resolution)
      # sum monthly precipitation (mm)
      # average temperatures (C)
      # average minimum temperatures
      # average maximum temperatures
      # max monthly snow water equivalent (mm)
 # ndvi rolling average 3, 6, 12, and 18 months
 # sum precipitation for 1, 3, 6, 12, 18 months previous
 # sum swe for 1, 3, 6, 12, 18 months previous
 # temp at time t, and rolling average of 3, 6, 12, and 18 months

 # test rolling average on a subset of data (zuni and navajo)
   pueblos.cov <- southwest %>%
      filter(site == "navajo" | site == "zuni")
      
    
 # rolling three-month average for ndvi
   pueblo.roll.ndvi <- pueblos.cov %>%
      group_by(site, web) %>%
      
      # rolling mean that month and the previous 3 months - worked
      mutate(ndvi3 = rollmean(ndvi, 4, align = "right", fill = NA)) %>%
      
      # rolling sum that month and the previous 3 months - worked
      mutate(prcp3  = rollsum(prcp, 4, align = "right", fill = NA))
   

 # worked - now apply to all 
   southwest.lags <- southwest %>%
      group_by(site, web) %>%
      
      # moving mean that month and the previous 3, 6, 12, and, 18months
      mutate(ndvi3    = rollmean(ndvi, 4,    align = "right", fill = NA)) %>%
      mutate(ndvi6    = rollmean(ndvi, 7,    align = "right", fill = NA)) %>%
      mutate(ndvi12   = rollmean(ndvi, 13,   align = "right", fill = NA)) %>%
      mutate(ndvi18   = rollmean(ndvi, 19,   align = "right", fill = NA)) %>%
      mutate(tmin3    = rollmean(tmin, 4,    align = "right", fill = NA)) %>%
      mutate(tmin6    = rollmean(tmin, 7,    align = "right", fill = NA)) %>%
      mutate(tmin12   = rollmean(tmin, 13,   align = "right", fill = NA)) %>%
      mutate(tmin18   = rollmean(tmin, 19,   align = "right", fill = NA)) %>%
      mutate(temp3    = rollmean(temp, 4,    align = "right", fill = NA)) %>%
      mutate(temp6    = rollmean(temp, 7,    align = "right", fill = NA)) %>%
      mutate(temp12   = rollmean(temp, 13,   align = "right", fill = NA)) %>%
      mutate(temp18   = rollmean(temp, 19,   align = "right", fill = NA)) %>%
      mutate(tmax3    = rollmean(tmax, 4,    align = "right", fill = NA)) %>%
      mutate(tmax6    = rollmean(tmax, 7,    align = "right", fill = NA)) %>%
      mutate(tmax12   = rollmean(tmax, 13,   align = "right", fill = NA)) %>%
      mutate(tmax18   = rollmean(tmax, 19,   align = "right", fill = NA)) %>%
      mutate(prcp3    = rollsum(prcp, 4,     align = "right", fill = NA)) %>%
      mutate(prcp6    = rollsum(prcp, 7,     align = "right", fill = NA)) %>%
      mutate(prcp12   = rollsum(prcp, 13,    align = "right", fill = NA)) %>%
      mutate(prcp18   = rollsum(prcp, 19,    align = "right", fill = NA)) %>%
      mutate(swe3     = rollsum(swe, 4,      align = "right", fill = NA)) %>%
      mutate(swe6     = rollsum(swe, 7,      align = "right", fill = NA)) %>%
      mutate(swe12    = rollsum(swe, 13,     align = "right", fill = NA)) %>%
      mutate(swe18    = rollsum(swe, 19,     align = "right", fill = NA)) 
   
   
 # UPDATED - 11/19/2019 - remove 18-month time lag (see folder README)
   updated.southwest.lags <- southwest %>%
      group_by(site, web) %>%
      
      # rolling mean that month and the previous 3, 6, and 12 months
      mutate(ndvi3    = rollmean(ndvi, 4,    align = "right", fill = NA)) %>%
      mutate(ndvi6    = rollmean(ndvi, 7,    align = "right", fill = NA)) %>%
      mutate(ndvi12   = rollmean(ndvi, 13,   align = "right", fill = NA)) %>%
      mutate(tmin3    = rollmean(tmin, 4,    align = "right", fill = NA)) %>%
      mutate(tmin6    = rollmean(tmin, 7,    align = "right", fill = NA)) %>%
      mutate(tmin12   = rollmean(tmin, 13,   align = "right", fill = NA)) %>%
      mutate(temp3    = rollmean(temp, 4,    align = "right", fill = NA)) %>%
      mutate(temp6    = rollmean(temp, 7,    align = "right", fill = NA)) %>%
      mutate(temp12   = rollmean(temp, 13,   align = "right", fill = NA)) %>%
      mutate(tmax3    = rollmean(tmax, 4,    align = "right", fill = NA)) %>%
      mutate(tmax6    = rollmean(tmax, 7,    align = "right", fill = NA)) %>%
      mutate(tmax12   = rollmean(tmax, 13,   align = "right", fill = NA)) %>%
      
      # sum over the past 3, 6, and 12 months
      mutate(prcp3    = rollsum(prcp, 4,     align = "right", fill = NA)) %>%
      mutate(prcp6    = rollsum(prcp, 7,     align = "right", fill = NA)) %>%
      mutate(prcp12   = rollsum(prcp, 13,    align = "right", fill = NA)) %>%
      mutate(swe3     = rollsum(swe, 4,      align = "right", fill = NA)) %>%
      mutate(swe6     = rollsum(swe, 7,      align = "right", fill = NA)) %>%
      mutate(swe12    = rollsum(swe, 13,     align = "right", fill = NA)) %>%
      
      # arrange
      select(date, month, site, web, 
             temp, temp3, temp6, temp12,
             tmin, tmin3, tmin6, tmin12,
             tmax, tmax3, tmax6, tmax12,
             prcp, prcp3, prcp6, prcp12,
             ndvi, ndvi3, ndvi6, ndvi12,
             swe,  swe3,  swe6,  swe12,
             swewinter)
      

   
## Correlation plot --------------------------------------------------------- ##
   
   
 # reorder dataframe columns
   southwest.order <- southwest.lags %>%
      dplyr::select(date, month, site, web,
                    tmin, tmin3, tmin6, tmin12, tmin18,
                    temp, temp3, temp6, temp12, temp18,
                    tmax, tmax3, tmax6, tmax12, tmax18,
                    prcp, prcp3, prcp6, prcp12, prcp18,
                    ndvi, ndvi3, ndvi6, ndvi12, ndvi18,
                    swe,  swe3,  swe6,  swe12,  swe18,
                    swewinter)
   
   
 # correlation plot of all the sites and covariates - subset to just covariates
   sw.subset.corr <- southwest.order %>%
      ungroup() %>%
      dplyr::select(-c(date, site, web))
   
   
 # format for a corrpot
   correlation <- cor(sw.subset.corr, use = "pairwise.complete.obs")
   
   
 # make a cor plot
   sw.correlation <- corrplot(correlation, method = "number")
    

## Center and scale data ---------------------------------------------------- ##
   
   
 # center and scale
   norm.southwest <- southwest.lags %>%
      mutate_at(scale, .vars = vars(5:35))
   
 # UPDATED - 11/19/2019 - normalize the new dataframe 
   updated.norm.southwest <- updated.southwest.lags %>%
      mutate_at(scale, .vars = vars(5:29))

   
 # write file
   write.csv(norm.southwest, "Data/southwest_covariates_norm.csv")
   
   
 # UPDATED - 11/19/2019 - save file in updated data folder
   write.csv(updated.norm.southwest, 
             "Data/Updated/updated_southwest_covariates_norm.csv")
   

 # UPDATED - 06/12/2020 - export unscaled covariates and combine dfs later
   write.csv(southwest.lags, 
             "Data/Updated/updated_southwest_covariates.csv")
   
   
 # save data
   save(southwest, 
        southwest.lags,
        norm.southwest,
        file = "Data/soutwest_covariates.RData")
   
 # UPDATED - 11/19/2019 - save covariates as RData
   save(southwest, 
        updated.southwest.lags,
        updated.norm.southwest,
        file = "Data/Updated/updated_soutwest_covariates.RData")

   
## Distributions of normalized data ----------------------------------------- ##
    

 # plot
   norm.southwest %>%
     keep(is.numeric) %>%
     gather() %>%
     ggplot(aes(value)) +
     facet_wrap(~ key, scales = "free") +
     geom_histogram()


## Correlation between variables -------------------------------------------- ##

   
 # update - print Pearson's cor coef
 

   
 # select variables - should i do this after center and scale?
   sw.merge.corr <- updated.norm.southwest %>%
      ungroup() %>%
     dplyr::select(-c(date, month, site, web))


 #  compute a correlation matrix
    cor.mat <- rcorr(as.matrix(sw.merge.corr, type = "pearson"))


 # make a correlogram - add month
   corrplot(cor.mat$r, method = "number", type = "upper", order = "hclust", 
            tl.col = "black", sig.level = 0.05, insig = "blank", tl.srt = 45)


 # chart correlation
   chart.Correlation(sw.merge.corr, histogram = TRUE, pch = 19)
     # distribution of each variable is shown on the diagonal
     # the bivariate scatter plots with a fitted line are displayed
     # the value of the correlation plus the significance level as stars


## Plot web and cov --------------------------------------------------------- ##


 # facet wrap for web and precipitation
   ggplot(southwest.order, aes(x = yearmon, y = prcp)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ web, ncol = 4)


 # facet wrap for web and ndvi
   ggplot(southwest.order, aes(x = yearmon, y = ndvi)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ web, ncol = 4)

 # facet wrap for web and swe
   ggplot(southwest.order, aes(x = yearmon, y = swe)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ web, ncol = 4)


 # facet wrap for web and tmax
   ggplot(ssouthwest.order, aes(x = yearmon, y = temp)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ web, ncol = 4)


## -------------------------------------------------------------------------- ##