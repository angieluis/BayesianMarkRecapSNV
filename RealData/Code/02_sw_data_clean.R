## Southwest data cleaning -------------------------------------------------- ##


 # source files
   source("Code/01_sw_data_functions.R")


## Read in combined capture data -------------------------------------------- ##


   southwest.dirty <- read.csv("Data/Covariates/southwest_data.csv", 
                               na.strings = c("", "NA"))
    # 8/3/2019 - changed from pema.dirty to southwest.dirty


## Organize ----------------------------------------------------------------- ##


 # keep pertinent columns
   pema <- southwest.dirty %>%
     dplyr::select(2:5, 9, 10, 13:17, 23:29, 32) %>%
     
     # mutate everything to lower case
     mutate_all(funs(tolower))
   
   
 # write to csv
   write.csv(pema, "Data/southwest_data_trim.csv")
    # 8/3/2019 - changed from pema.csv to southwest.csv
   
   
 # change some columns to factors and format date
   pema$tag      <- as.factor(pema$tag)
   pema$letter_2 <- as.character(pema$letter_2)
   pema$sex      <- as.factor(pema$sex)
   pema$snv_pos  <- as.factor(pema$snv_pos)
   pema$date     <- mdy(pema$date)          


 # make sure to remove the nine NA rows at the bottom of the dataframe
   pema <-pema[-c(84444, 84445, 84446, 84447, 84448, 84449, 84450, 84451, 84452), ]    


## Make dataframe with lags times between each date ------------------------- ##


 # make a new column that adds unique id to each session (numeric)
   pema.lag <- pema %>%
     
     # paste site and session together
     mutate(site_session = paste(site, session, sep = "_")) %>%
  
     # create an id (starting at 1) for each session based on site
     group_by(site) %>%
  
     # omg this took altogether too long
     mutate(session_id = match(site_session, unique(site_session)))


   # change date from factor
     #pema.lag$date <- as.POSIXct(pema.lag$date, format = "%Y-%m-%d")


 # make new variable lag.date 
   lag.date <- pema.lag %>%
     
     # get distinct dates and session_id
     distinct(site, date, session_id) %>%
  
     # sort by site then date
     arrange(site, date, session_id) %>%
  
     # select date and session_id - also automatically selects site...
     dplyr::select(site, date, session_id) %>%
  
     # make new column, diff_date - first day will be NA since no date prior
     mutate(diff_date = c(NA, diff(date))) %>%
  
     # round
     mutate_if(is.numeric, round, 1)


## Drop rows ---------------------------------------------------------------- ##


 pema.pre.dupe <- pema.lag %>%
   
   # column to tie tags to site (lots of tag # repeats across sites)
   mutate(tag_site = paste(tag, site, sep = "_")) %>%
  
   # column to tie webs to sites (some webs have the same name)
   mutate(site_web = paste(site, web, sep = "_")) %>%

   # drop dead in trap/on table, -9 tag, and refugia
   filter(!fate %in% c(4, 14, 24, 34), 
          !tag == -9, 
          !site == "refugia",
          !site == "sevilleta") %>%
   
   # change juveniles that test positive to negative (maternal antibodies)
   mutate(snv_pos = replace(snv_pos, age == 1 & snv_pos == 1, 0)) %>%

   # make one reproductive column - yes/no
   mutate(repro = if_else(testes == 0 |
                            vagina == 0 |
                            nips_en == 0 |
                            nips_ln == 0 |
                            pregnant == 0, 0, if_else(testes == 1 |
                                    vagina == 1 |
                                    nips_en == 1 |
                                    nips_ln == 1 |
                                    pregnant == 1, 1, -9))) %>%
   
   # make one wound column - yes/no
   mutate(wound = if_else(wound == 1 |
                            sever == 1 |
                            sever == 2 |
                            sever == 3, 1,
                          if_else(wound == 0 |
                                    sever == 0, 0, -9))) %>%

   # drop old repro columns and wound columns
   dplyr::select(-c(testes, vagina, nips_en, nips_ln, pregnant, sever)) %>%
  
   # rearrange columns...sigh...drop age, repro, wound, and mass
   dplyr::select(univ, site, web, site_web, session, session_id, site_session, date,
          tag, tag_site, fate, letter_2, sex, snv_pos) %>%
  
   # mutate everything to lowercase
   mutate_all(funs(tolower)) 


 # change date from character to date
   pema.pre.dupe$date     <- as.Date(pema.pre.dupe$date, format = "%Y-%m-%d") 


## Look for duplicates ------------------------------------------------------ ##


 # made a function that evaluates duplicates see dupe.fun
 # some animals showing up same tag, same day
   dupes <- dupes.fun(pema.pre.dupe, date, tag, letter_2)


 # arrange by date
   dupes <- dupes %>%
     
     # group by site
     group_by(site) %>%
     
     # the arrange in descending order
     arrange(desc(date), .by_group = TRUE)
   

 # ~44 instances of duplicates (i.e., same tag same day), differences range from
 # being on different webs to having different sex. function does not show first 
 # occurance


 # find duplicates and decide what case to remove based on first occurrance
 # add values based on pema.dupes - ugh, case by case


 # dupes with no dates but contain different values for certain variables
 # remove site == "molina", tag == "1262",  web == "mb", session_id == "17"
 # remove site == "navajo", tag == "6604",  web == "2",  session_id == "50"
 # remove site == "zuni", tag == "15215", snv_pos == "-9, session_id == "123"
 # remove site == "navajo", tag == "3239",  web == "1"
 # remove site == "grandcanyon", tag == "i931",  web == "e"
 # remove site == "grandcanyon", tag == "i640",  web == "e"
 # remove site == "hesperus", tag == "1706",  web == "ha"
 # remove site == "placitas", tag == "10932", web == "3"
 # remove site == "zuni", tag == "15665,  web == "3"
 # remove site == "zuni", tag == "15427", web == "3"
 # remove site == "hesperus", tag == "2503", sex == "0"
 # remove site == "placitas", tag == "tc040(25)p", web == "3"
 # remove site == "placitas", tag == "tc(32)500", web == "3" 

 # dupes with dates involved
 # remove site == "placitas", tag == "10812", fate == "3", date == "3/22/2005"
 # remove site == "zuni", tag == "15266", fate == "1", date == "10/12/2005"
 # remove site == "zuni", tag == "15897", fate == "3", date == "1/10/2006"

 # dupe sites with no differences
 # site == "placitas", tag == "10887", date == "3/24/2005"
 # site == "navajo",   tag == "14034", date == "8/11/2005"
 # site == "zuni",     tag == "15224", date == "10/11/2005"
   
 # step one - remove weird ass duplicates
 pema.adj <- pema.pre.dupe %>% 
   
   # filter on specific conditions without dates
   filter(!site == "navajo" | !tag == "3239" | !web == "1") %>%
   filter(!site == "grandcanyon" | !tag == "i931" | !web == "e") %>%
   filter(!site == "grandcanyon" | !tag == "i640" | !web == "e") %>%
   filter(!site == "hesperus" | !tag == "1706" | !web == "ha") %>%
   filter(!site == "placitas" | !tag == "10932" | !web == "3") %>%
   filter(!site == "zuni" | !tag == "15427" | !web == "3") %>%
   filter(!site == "hesperus" | !tag == "2503" | !sex == "0") %>%
   filter(!site == "placitas" | !tag == "tc040(25)p" | !web == "3") %>%
   filter(!site == "placitas" | !tag == "tc(32)500" | !web == "3") %>%
   filter(!site == "zuni" | !tag == "15665" | !web == "3") %>%
  
   # step two - filter dupes where session is a criteria
   filter(!site == "molina" | !tag == "1262" | !web == "mb" | !session_id == "17") %>%
   filter(!site == "navajo" | !tag == "6604" | !web == "2" | !session_id == "50") %>%
   filter(!site == "navajo" | !tag == "6604" | !web == "2" | !session_id == "50") %>%
   filter(!site == "zuni" | !tag == "15215" | !snv_pos == "-9" | !session_id == "123") %>%
  
   # step three - filter dupes where date is a criteria
   filter(!site == "placitas" | !tag == "10812" | !fate == "3" | !date == "2005-03-22") %>%
   filter(!site == "zuni" | !tag == "15266" | !fate == "1" | !date == "2005-10-12") %>%
   filter(!site == "zuni" | !tag == "15897" | !fate == "3" | !date == "2006-01-10") %>%
   filter(!site == "pcms" | !tag == "1846" | !fate == "1" | !date == "1997-07-15")


 # step four -  remove the true duplicates (i.e., everything same)
   pema.adj <- distinct(pema.adj)
   # yaaaas


## Handle NA and -9 data for letter_2 --------------------------------------- ##


 # a lot of -9 and NA - lots of info not input if it was a recapture 
 ind.sp <- which(pema.adj$letter_2 == -9 |
                   is.na(pema.adj$letter_2) |
                   pema.adj$letter_2 == "")


 # test ind for the mentioned logical contraints
   recap <- logical()


 # loop through individuals to see if -9, NA, or "" were recaps
   for (i in 1:length(ind.sp)) {
     recap[i] <- length(which(pema.adj$tag_site == 
                                pema.adj$tag_site[ind.sp[i]])) > 1
     }


 # delete if no species information (w/ no recaps to look under)
   pema.adj <- pema.adj[-ind.sp[which(recap == FALSE)], ]


 # look for species info when there are recaps for the tag number
   ind.sp <- which(pema.adj$letter_2 == -9 |
                     is.na(pema.adj$letter_2) |
                     pema.adj$letter_2 == "")


 # test character type objects
   sp <- character()


 # loop through tag_site
   for (i in 1:length(ind.sp)) { 
     sp1   = as.character(pema.adj$letter_2[which(pema.adj$tag_site ==
                                                     pema.adj$tag_site[ind.sp[i]])])
      # return the first thing that isn't NA or -9
      sp[i] = sp1[!is.na(sp1) & sp1 != -9][1]
      pema.adj$letter_2[ind.sp[i]] = sp[i]
      }


 # remove recaps that had no species id
   ind.sp <- which(pema.adj$letter_2 == -9 |
                     is.na(pema.adj$letter_2) | 
                     pema.adj$letter_2 == "")


 # remove the last individuals with no information 
   pema.adj <- pema.adj[-ind.sp,]


 # NOTE - worked
   which(summary(factor(pema.adj$letter_2)) > 3000)
   unique(pema.adj$letter_2)


## Fill sex based on previous captures -------------------------------------- ##


 # a lot of -9 and NA, lots of info wasn't input if it was a recapture (4281)
   ind.sex <- which(pema.adj$sex == -9 | 
                      is.na(pema.adj$sex) | 
                      pema.adj$sex== "")


 # test ind for the mentioned logical contraints
   recap <- logical()


 # loop through individuals to see if -9, NA, or "" were recaps
   for (i in 1:length(ind.sex)) {
     recap[i] = length(which(pema.adj$tag_site == 
                               pema.adj$tag_site[ind.sex[i]])) > 1
     }


 # delete if no species information (w/ no recaps to look under)
   pema.adj <- pema.adj[-ind.sex[which(recap == FALSE)], ]


 # look for species info when there are recaps for the tag number
   ind.sex <- which(pema.adj$sex == -9 |
                      is.na(pema.adj$sex) |
                      pema.adj$sex == "")


 # test character type objects
   sp <- character()


 # loop through tag_site
   for (i in 1:length(ind.sex)) {
     sp1   = as.character(pema.adj$sex[which(pema.adj$tag_site ==
                                               pema.adj$tag_site[ind.sex[i]])])
     # return the first thing that isn't NA or -9
     sp[i] = sp1[!is.na(sp1) & sp1 != -9][1]
     pema.adj$sex[ind.sex[i]] = sp[i]
     } # this loop backfilled all the data for -9, left ~68 NAs


 # NOTE - worked but still 26 NA's
   summary(factor(pema.adj$sex))


## Fill in serostatus of adults based on previous captures ------------------ ##


 # how many -9 in pema.adj$snv_pos before we run loop?
 # -2 = 5, -3 = 3050, -9 = 1753, 0 = 47779, 1 = 5420, NA = 4656
   summary(factor(pema.adj$snv_pos))


 # find lines that have missing snv data
   ind.snv <- which(pema.adj$snv_pos == -9 |
                      pema.adj$snv_pos == -3 |
                      pema.adj$snv_pos == "" |
                      is.na(pema.adj$snv_pos))


 # create an object for length of snv
   l <- numeric()


 # for loop to move through individuals by tag_site and session
   for (i in 1:length(ind.snv)) {
     tag_site = pema.adj$tag_site[ind.snv[i]]
     session  = pema.adj$session[ind.snv[i]]
     snv      = pema.adj$snv_pos[which((pema.adj$tag_site == tag_site &
                                        pema.adj$session == session))]
     l[i]     = length(snv)
     # get all snv statuses for that individual, that month
     # if there is data fill in values
     if (length(snv) > 1) {
       pema.adj$snv_pos[ind.snv[i]] = snv[which(snv == 1 | snv == 0)][1]
     }
     }


 # NOTE - worked but added ~4600 NA's
   summary(factor(pema.adj$snv_pos))


## Change all values but 1 or 0 to NA --------------------------------------- ##


 pema.est <-  pema.adj %>%
   
   # replace values that aren't 0 or 1 to NA - to estimate infection later
   # for the non-MS model change NAs to different value so as to not confuse
   # the model (NAs not trapped)
   # 3 means unknown serostatus
   mutate(snv_adj = replace(snv_pos, which(snv_pos == -9 |
                                             snv_pos == -3 | 
                                             snv_pos == -2 |
                                             is.na(snv_pos)), NA))


 # NOTE - worked ~9464 with unknown infection status
   summary(factor(pema.est$snv_adj))
   
   
 # how many positive males - 671 males
   male.pos    <- pema.est %>%
     distinct(tag_site, .keep_all = TRUE) %>%
     filter(sex == 1 & snv_pos == 1)
   
   
 # how many positive females - 212 females
   fem.pos    <- pema.est %>%
     distinct(tag_site, .keep_all = TRUE) %>%
     filter(sex == 0 & snv_pos == 1)
   


## Look for duplicates again ------------------------------------------------ ##

 
   dupes <- dupes.fun(pema.est, date, site, tag, letter_2)
    # worked - no duplicates


## Look to see if an animal moved between sites ----------------------------- ##


   move <- pema.est %>%
     dplyr::select(tag_site, web) %>%
     group_by(tag_site) %>%
     summarise(n_distinct(web))
  

 # 307 out of 21331 occurences where an animal "used" a different web
   sum(move$`n_distinct(web)` > 1)
   sum(move$`n_distinct(web)`)


## Were webs trapped on the same date --------------------------------------- ##


 # group by site and session to get first trapping date within a session
   first.date <- pema.est %>%
     dplyr::select(site_web, session, date) %>%
     group_by(site_web, session) %>%
     distinct(date) %>%
     arrange(date) %>%
     slice(1) %>%
     spread(site_web, date)


 # write table
   write.table(first.date, file = "Data/southwest_first_date.csv", sep = ",")
   
   
## -------------------------------------------------------------------------- ##
   
   
 # final cleaned data name - all species and all sites
   southwest.final.clean <- pema.est
    # 8/3/2019 - changed from pema.final.clean to southwest.final.clean
   

## -------------------------------------------------------------------------- ##
 
   
 # save dirty and clean dataframes 
   save(southwest.dirty, 
        southwest.final.clean, 
        file = "Data/AllCaptureData.RData")
   

## -------------------------------------------------------------------------- ##