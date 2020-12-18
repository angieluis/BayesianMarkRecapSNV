# Navajo - Web 1 & Web 2: Visualization ------------------------------------- ##
## -------------------------------------------------------------------------- ##


  # NOTES
  # variables included in Navajo model
   # ndvi0     (ndvi average at a 0-month lag)
   # tmin0     (temperature min at a 0-month lag)
   # tmax0     (temperature max at a 0-month lag)
   # swe       (snow-water equiv. at a 0-month lag)
   # swewinter (snow-water equiv. sum for spring/summer - vector of 7)
    # oct-mar
    # oct-apr
    # oct-may
    # oct-jun
    # oct-jul
    # oct-aug
    # oct-sep


 # packages
   library(coda)
   library(MCMCvis)
   library(mcmcplots)
   library(tidyverse)
   library(ggmcmc)
   library(lubridate)
   library(reshape2)
   library(ggridges)
   library(ggeffects)
   library(viridis)
   library(scales)

   
 # load RData specific to site
   load("Data/ModelRuns/ForAnalysis/Navajo_MaxModel_June2020.RData")
   
   
 # read in scaled covariate data
   scale_covariates <- read.csv("Data/Updated/updated_southwest_covariates_norm.csv")
   
   
 # read in unscaled covariate data
   unscale_covariates <- read.csv("Data/Updated/updated_southwest_covariates.csv")
   
   
 # navajo
   navajo_model <- N12.rCJS.MaxModel.J22
   # navajo_model_update <- update.N12.rCJS.MaxModel
   navajo_sims  <- N12.rCJS.MaxModel.J22$BUGSoutput$sims.list 
   
   
 # specify parameters
   params <-  c("alpha.0",
                "alpha.male", 
                "alpha.web", 
                "alpha.swe",
                "alpha.swe.winter",
                "alpha.season",
                "alpha.ndvi.season",
                "alpha.tmin.season",
                "alpha.tmax.season", 
                "sigma.0", 
                "sigma.male",
                "sigma.recap",
                "sigma.season")
   
   
## Functions for visualization ---------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
 # function that takes model of interest and displays various outputs
   bugs_output_fun <- function(model) {
      
      
      # require MCMCvis for visualization
      require(MCMCvis)
      require(mcmcplots)
      require(tidyverse)
      
      
      # model summary
      model.summary = MCMCsummary(model,
                                  round = 4)
      
      
      # plot detailed
      # revisit parameters needed as you run through models
      model.plot.detail = MCMCplot(object  = model,
                                   params  = params,
                                   xlim    = c(-4, 4),
                                   ref_ovl = TRUE,
                                   ISB = FALSE)
      
      # plot mcmcplots output
      mcmcplot = mcmcplot(model,
                          parms = params)
    
      # output of model
      output = list(model.summary,
                    model.plot.detail,
                    mcmcplot
                    
      )
      
      
      # returns the summary and ind.results in console
      # trace and plots pop up automatically
      return(output)
      
   } 

   
 # look at basic model output function
   bugs_output_fun(navajo_model) 


 # reverse logit function
   rev.logit <- function(x) {
      exp(x) / (1 + exp(x))
   }

   
## UPDATE: MCMCplot for thesis ---------------------------------------------- ##
  
   
 # summary and posterior estimates
   navajo_summary <- MCMCsummary(navajo_model, params = params)
 
   
 # write both to a csv
   #write.csv(navajo_summary, "navajo_summary.csv")

   

 # calculate prior posterior overlap
   PR <- rnorm(15000, 0, 32) # equivalent to dnorm(0, 0.001) in JAGS
   PPO <- MCMCtrace(navajo_model, 
                    params  = params,
                    ISB     = FALSE,
                    priors  = PR,
                    plot    = FALSE,
                    PPO_out = TRUE)
   
   
 # update posterior plot
   navajo_mcmc <- MCMCplot(object      = navajo_model,
                           params      = params,
                           labels      = c("phi - intercept",
                                           "sex", 
                                           "web", 
                                           "swe",
                                           "swewinter",
                                           "spring",
                                           "summer",
                                           "autumn",
                                           "ndvi0 spring",
                                           "ndvi0 summer",
                                           "ndvi0 autumn",
                                           "tmin0 spring",
                                           "tmin0 summer",
                                           "tmin0 autumn",
                                           "tmax0 spring",
                                           "tmax0 summer",
                                           "tmax0 autumn",
                                           "p - intercept", 
                                           "sex",
                                           "recapture",
                                           "spring",
                                           "summer",
                                           "autumn"),
                           rank        = FALSE, 
                           xlim        = c(-4, 4),
                           xlab        = "Navajo - Estimates",
                           guide_lines = TRUE,
                           ref_ovl     = TRUE,
                           ISB         = FALSE)
   
   
   
   
## Plot Navajo phi - these are relatively useless ----------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
   
 # navajo - female web 1
   navajo_f_w1 <- navajo_sims$phi.female.web1 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nf1" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - female web 2
   navajo_f_w2 <- navajo_sims$phi.female.web2 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nf2" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - male web 1
   navajo_m_w1 <-navajo_sims$phi.male.web1 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nm1" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - male web 2
   navajo_m_w2 <- navajo_sims$phi.male.web2 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nm2" = ".") %>%
      rowid_to_column("date.id")
   
   
 # join dataframs
   navajo_phi_merge <- navajo_f_w1 %>%
      left_join(navajo_f_w2) %>%
      left_join(navajo_m_w1) %>%
      left_join(navajo_m_w2)
   
   
 # melt data
   navajo_phi_melt <- melt(navajo_phi_merge, id.vars = "date.id")
   
   
 # plot
   navajo_phi_plot <- ggplot(navajo_phi_melt,
                             aes(x     = variable, 
                                 y     = value,
                                 fill  = variable)) +
      
      stat_boxplot(geom = "errorbar") +
      
      geom_boxplot() +
      
      scale_fill_viridis(discrete = TRUE) +
      
      geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_discrete(name   = "\nGroup",
                       labels = c("Female Web1",
                                  "Female Web2",
                                  "Male Web1",
                                  "Male Web2")) +
      
      ggtitle("Mean Survival - Navajo\n") +
      
      theme_bw() +
      
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black")) 
 
   
## Plot Navajo p -------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##  
   
   
 # navajo - female web 1
   navajo_fp_w1 <- navajo_sims$p.female.web1 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nf1.p" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - female web 1
   navajo_fp_w2 <- navajo_sims$p.female.web2 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nf2.p" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - male web 1
   navajo_mp_w1 <-navajo_sims$p.male.web1 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nm1.p" = ".") %>%
      rowid_to_column("date.id")
   
   
 # navajo - male web 2
   navajo_mp_w2 <- navajo_sims$p.male.web2 %>% 
      colMeans() %>%
      as.data.frame() %>%
      rename("nm2.p" = ".") %>%
      rowid_to_column("date.id")
   
   
 # join dataframs
   navajo_p_merge <- navajo_fp_w1 %>%
      left_join(navajo_fp_w2) %>%
      left_join(navajo_mp_w1) %>%
      left_join(navajo_mp_w2)
   
   
 # melt data
   navajo_p_melt <- melt(navajo_p_merge, id.vars = "date.id")
   
   
 # plot
   navajo_p_plot <- ggplot(navajo_p_melt, 
                           aes(x    = variable, 
                               y    = value,
                               fill = variable)) + 
      
      stat_boxplot(geom = "errorbar") +
      
      geom_boxplot() +
      
      geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
      
      scale_fill_viridis(discrete = TRUE) +
      
      scale_y_continuous(name   = "Recapture\n",
                         breaks = seq(0.2, 0.25, 0.05),
                         limits = c(0.2, 0.25)) +
      
      scale_x_discrete(name   = "\nGroup",
                       labels = c("Female Web1",
                                  "Female Web2",
                                  "Male Web1",
                                  "Male Web2")) +
      
      ggtitle("Mean Recapture - Navajo\n") +
      
      theme_bw() +
      
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
   

## Calculate means and spread for environmental variables ------------------- ##
## -------------------------------------------------------------------------- ##
   
   
 # group by site_web and take and season (vector of 4)
    # 1 = december - march: accounts for ENSO events (12, 1, 2, 3)
    # 2 = april - june (4, 5, 6)
    # 3 = july - september: accounts for monsoon events (7, 8, 9)
    # 4 = october - november (10, 11)
   
 # swewinter (snow-water equiv. sum for spring/summer - vector of 7)
    # oct-mar - season 1
    # oct-apr - season 2
    # oct-may - season 2
    # oct-jun - season 2
    # oct-jul - season 3
    # oct-aug - season 3
    # oct-sep - season 3
 
   
## Combine scaled and unscaled data ----------------------------------------- ##
   
  
 # check dims of dataframes
   
 # 3528x36
   dim(unscale_covariates)
   
 # 3528x30
   dim(scale_covariates)
   
   
 # drop unnecessary columns and rename unscaled data
   unscaled <- unscale_covariates %>%
      filter(site == "navajo") %>%
      filter(web == "1") %>%
      select(date, ndvi, tmin, tmax, swe, swewinter) %>%
      rename(ndvi_raw      = ndvi,
             tmin_raw      = tmin,
             tmax_raw      = tmax,
             swe_raw       = swe,
             swewinter_raw = swewinter)
   
   
 # same ish with scaled covariates
   scaled <- scale_covariates %>%
      filter(site == "navajo") %>%
      filter(web == "1") %>%
      select(date, month, site, web, ndvi, tmin, tmax, swe, swewinter)
   
   
 # combine the dataframes
   cov_joined <- scaled %>%
      left_join(unscaled, by = "date")
   

## Web 1 covariate data ----------------------------------------------------- ##
   
   
 # subset web 1
   cov_web1 <- scale_covariates %>%
      filter(site == "navajo") %>%
      filter(web == "1") %>%
      select(date, month, site, web, ndvi, tmin, tmax, swe, swewinter)
   
   
 # group winter and take the mean of variables
   winter_web1 <- cov_web1 %>%
      filter(month == 1 | month == 2 | month == 3 | month == 12) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group spring and take the mean of variables
   spring_web1 <- cov_web1 %>%
      filter(month == 4 | month == 5 | month == 6) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group summer and take the mean of variables
   summer_web1 <- cov_web1 %>%
      filter(month == 7 | month == 8 | month == 9) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group winter and take the mean of variables
   autumn_web1 <- cov_web1 %>%
      filter(month == 10 | month == 11) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # values for mean 
   mean.ndvi.w1       <- c(winter_web1$ndvi, 
                           spring_web1$ndvi,
                           summer_web1$ndvi,
                           autumn_web1$ndvi)
   
   
   mean.tmin.w1       <- c(winter_web1$tmin, 
                           spring_web1$tmin,
                           summer_web1$tmin,
                           autumn_web1$tmin)
   
   
   mean.tmax.w1       <- c(winter_web1$tmax, 
                           spring_web1$tmax,
                           summer_web1$tmax,
                           autumn_web1$tmax)
   
   
   mean.swe.w1        <-  c(winter_web1$swe, 
                            spring_web1$swe,
                            summer_web1$swe,
                            autumn_web1$swe)
   
   
   mean.swe.winter.w1 <-  c(winter_web1$swewinter, 
                            spring_web1$swewinter,
                            summer_web1$swewinter,
                            autumn_web1$swewinter)
   
   
 # values for spread calculations
   winter_ndvi_web1_sp <- winter_web1$ndvi
   spring_ndvi_web1_sp <- spring_web1$ndvi
   summer_ndvi_web1_sp <- summer_web1$ndvi
   autumn_ndvi_web1_sp <- autumn_web1$ndvi
   
   
   winter_tmin_web1_sp <- winter_web1$tmin
   spring_tmin_web1_sp <- spring_web1$tmin
   summer_tmin_web1_sp <- summer_web1$tmin
   autumn_tmin_web1_sp <- autumn_web1$tmin
   
   
   winter_tmax_web1_sp <- winter_web1$tmax
   spring_tmax_web1_sp <- spring_web1$tmax
   summer_tmax_web1_sp <- summer_web1$tmax
   autumn_tmax_web1_sp <- autumn_web1$tmax
   
   
   winter_swe_web1_sp <- winter_web1$swe 
   spring_swe_web1_sp <- spring_web1$swe
   summer_swe_web1_sp <- summer_web1$swe
   autumn_swe_web1_sp <- autumn_web1$swe
   
   
   winter_swewinter_web1_sp <- winter_web1$swewinter
   spring_swewinter_web1_sp <- spring_web1$swewinter
   summer_swewinter_web1_sp <- summer_web1$swewinter
   autumn_swewinter_web1_sp <- autumn_web1$swewinter
   
   
## Web 2 covariate data ----------------------------------------------------- ##
   
   
 # subset web 2 - use scaled
   cov_web2 <- scale_covariates %>%
      filter(site == "navajo") %>%
      filter(web == "2") %>%
      select(date, month, site, web, ndvi, tmin, tmax, swe, swewinter)
   
   
 # group winter and take the mean of variables
   winter_web2 <- cov_web2 %>%
      filter(month == 1 | month == 2 | month == 3 | month == 12) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group spring and take the mean of variables
   spring_web2 <- cov_web2 %>%
      filter(month == 4 | month == 5 | month == 6) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group summer and take the mean of variables
   summer_web2 <- cov_web2 %>%
      filter(month == 7 | month == 8 | month == 9) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # group winter and take the mean of variables
   autumn_web2 <- cov_web2 %>%
      filter(month == 10 | month == 11) %>%
      summarize_at(c("ndvi", "tmin", "tmax", "swe", "swewinter"), 
                   mean, na.rm = TRUE)
   
   
 # assign values for each season
   mean.ndvi.w2       <- c(winter_web2$ndvi, 
                           spring_web2$ndvi,
                           summer_web2$ndvi,
                           autumn_web2$ndvi)
   
   
   mean.tmin.w2       <- c(winter_web2$tmin, 
                           spring_web2$tmin,
                           summer_web2$tmin,
                           autumn_web2$tmin)
   
   
   mean.tmax.w2       <- c(winter_web2$tmax, 
                           spring_web2$tmax,
                           summer_web2$tmax,
                           autumn_web2$tmax)
   
   
   mean.swe.w2        <-  c(winter_web2$swe, 
                            spring_web2$swe,
                            summer_web2$swe,
                            autumn_web2$swe)
   
   
   mean.swe.winter.w2 <-  c(winter_web2$swewinter, 
                            spring_web2$swewinter,
                            summer_web2$swewinter,
                            autumn_web2$swewinter)
   
   
 # vectorize values for spread
   winter_ndvi_web2_sp <- winter_web2$ndvi
   spring_ndvi_web2_sp <- spring_web2$ndvi
   summer_ndvi_web2_sp <- summer_web2$ndvi
   autumn_ndvi_web2_sp <- autumn_web2$ndvi
   
   
   winter_tmin_web2_sp <- winter_web2$tmin
   spring_tmin_web2_sp <- spring_web2$tmin
   summer_tmin_web2_sp <- summer_web2$tmin
   autumn_tmin_web2_sp <- autumn_web2$tmin
   
   
   winter_tmax_web2_sp <- winter_web2$tmax
   spring_tmax_web2_sp <- spring_web2$tmax
   summer_tmax_web2_sp <- summer_web2$tmax
   autumn_tmax_web2_sp <- autumn_web2$tmax
   
   
   winter_swe_web2_sp <- winter_web2$swe 
   spring_swe_web2_sp <- spring_web2$swe
   summer_swe_web2_sp <- summer_web2$swe
   autumn_swe_web2_sp <- autumn_web2$swe
   
   
   winter_swewinter_web2_sp <- winter_web2$swewinter
   spring_swewinter_web2_sp <- spring_web2$swewinter
   summer_swewinter_web2_sp <- summer_web2$swewinter
   autumn_swewinter_web2_sp <- autumn_web2$swewinter
   
   
## Set values for mean estimates -------------------------------------------- ##
## -------------------------------------------------------------------------- ## 
   
   
 # use the mean of the values - look at spread later
   alpha.0          <- mean(navajo_sims$alpha.0)
   alpha.male       <- mean(navajo_sims$alpha.male)
   alpha.season.use <- c(0, 
                         mean(navajo_sims$alpha.season[ , 1]),
                         mean(navajo_sims$alpha.season[ , 2]),
                         mean(navajo_sims$alpha.season[ , 3]))
   
   alpha.web.use    <- c(0, 
                         navajo_sims$alpha.web[ , 1])
   
   
 # coefficients for environmental covariates
  # alpha.ndvi       <- mean(navajo_sims$alpha.ndvi)
  # alpha.tmin       <- mean(navajo_sims$alpha.tmin)
  # alpha.tmax       <- mean(navajo_sims$alpha.tmax)
   alpha.swe        <- mean(navajo_sims$alpha.swe)
   alpha.swe.winter <- mean(navajo_sims$alpha.swe.winter)
   
   
 # seasonal interactions
   alpha.ndvi.season.use <- c(0, 
                              mean(navajo_sims$alpha.ndvi.season[ , 1]),
                              mean(navajo_sims$alpha.ndvi.season[ , 2]),
                              mean(navajo_sims$alpha.ndvi.season[ , 3]))
   
   
   alpha.tmin.season.use <- c(0, 
                              mean(navajo_sims$alpha.tmin.season[ , 1]),
                              mean(navajo_sims$alpha.tmin.season[ , 2]),
                              mean(navajo_sims$alpha.tmin.season[ , 3]))
   
   
   alpha.tmax.season.use <- c(0, 
                              mean(navajo_sims$alpha.tmax.season[ , 1]),
                              mean(navajo_sims$alpha.tmax.season[ , 2]),
                              mean(navajo_sims$alpha.tmax.season[ , 3]))

   
## Calculate mean seasonal dynamics  ---------------------------------------- ##
## -------------------------------------------------------------------------- ##
 
   
 # QUESTION
 # We didn't set up the equation correctly in this context and I keep
 # receiving this error:
 # "longer object length is not a multiple of shorter object length"
   # Check - alpha.season.use or alpha.web.use (I think that's the issue)

   
 # specify equations - should return a vector of length 4, one for each season
 # equation for female web 1
   phi.female.web1.seasonal <- rev.logit(
      alpha.0 +
         alpha.male * 0 +
         alpha.season.use + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1 +
         # alpha.tmin * mean.tmin.w1 + 
         # alpha.tmax * mean.tmax.w1 +
         alpha.swe * mean.swe.w1 + 
         alpha.swe.winter * mean.swe.winter.w1 +
         alpha.ndvi.season.use * mean.ndvi.w1 +
         alpha.tmin.season.use * mean.tmin.w1 +
         alpha.tmax.season.use * mean.tmax.w1
   )
   
   
 # equation for female web 2
   phi.female.web2.seasonal <- rev.logit(
      alpha.0 +
         alpha.male * 0 +
         alpha.season.use + 
         alpha.web.use[2] + 
         # alpha.ndvi * mean.ndvi.w2 +
         # alpha.tmin * mean.tmin.w2 + 
         # alpha.tmax * mean.tmax.w2 +
         alpha.swe * mean.swe.w2 + 
         alpha.swe.winter * mean.swe.winter.w2 + 
         alpha.ndvi.season.use * mean.ndvi.w2 +
         alpha.tmin.season.use * mean.tmin.w2 +
         alpha.tmax.season.use * mean.tmax.w2
   )
   
   
 # equation for male web 1
   phi.male.web1.seasonal <- rev.logit(
      alpha.0 + 
         alpha.male * 1 + 
         alpha.season.use + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1 +
         # alpha.tmin * mean.tmin.w1 + 
         # alpha.tmax * mean.tmax.w1 + 
         alpha.swe * mean.swe.w1 + 
         alpha.swe.winter * mean.swe.winter.w1 + 
         alpha.ndvi.season.use * mean.ndvi.w1 +
         alpha.tmin.season.use * mean.tmin.w1 +
         alpha.tmax.season.use * mean.tmax.w1
   )
   
   
 # equation for male web 2
   phi.male.web2.seasonal <- rev.logit(
      alpha.0 + 
         alpha.male * 1 + 
         alpha.season.use + 
         alpha.web.use[2] + 
         # alpha.ndvi * mean.ndvi.w2 +
         # alpha.tmin * mean.tmin.w2 + 
         # alpha.tmax * mean.tmax.w2 + 
         alpha.swe * mean.swe.w2 + 
         alpha.swe.winter * mean.swe.winter.w2 + 
         alpha.ndvi.season.use * mean.ndvi.w2 +
         alpha.tmin.season.use * mean.tmin.w2 +
         alpha.tmax.season.use * mean.tmax.w2
   )

   
## Plot seasonal mean (with error??) ---------------------------------------- ##
   
   
 # make numeric output into dataframe
   fw1.seasonal <- as.data.frame(phi.female.web1.seasonal)
   mw1.seasonal <- as.data.frame(phi.male.web1.seasonal)
   fw2.seasonal <- as.data.frame(phi.female.web2.seasonal)
   mw2.seasonal <- as.data.frame(phi.male.web2.seasonal)
   
   
 # bind dataframes
   mean.season <- cbind(fw1.seasonal,
                        mw1.seasonal,
                        fw2.seasonal,
                        mw2.seasonal)
   
   
 # rename columns
   mean.season <- mean.season %>%
      rename(fw1.seasonal = phi.female.web1.seasonal,
             mw1.seasonal = phi.male.web1.seasonal,
             fw2.seasonal = phi.female.web2.seasonal,
             mw2.seasonal = phi.male.web2.seasonal) %>%
      
      # make a date id column
      rowid_to_column("date.id")
   
   
 # melt data frame ()
   mean.season.melt <- melt(mean.season, id.vars = "date.id")
   
   
 # plot a connected -  try lollipop 
   mean.season.plot <- ggplot(mean.season.melt,
                              aes(x    = as.factor(date.id),
                                  y    = value,
                                  color = variable,
                                  group = variable)) +
      
      geom_line() + 
      
      geom_point() +
      
      # geom_bar(stat     = "identity", 
       #        width    = 0.5, 
        #       position = "dodge",
         #      color    = "black") +
    
      scale_color_viridis(name     = "Group",
                          labels   = c("Female Web1",
                                       "Male Web1",
                                       "Female Web2",
                                       "Male Web2"),
                          discrete = TRUE) +
      
      scale_x_discrete(name   = "\nSeason",
                       breaks = c("1", "2", "3", "4"),
                       labels = c("Winter",
                                  "Spring",
                                  "Summer",
                                  "Autumn")) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0.0, 1.0, 0.1),
                         limits = c(0.0, 1.0)) +
      
      ggtitle("Mean Survival by Season - Navajo\n") +
      
      theme_bw() +
      
      theme(panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
   
   

## Set values for spread (violin plot) -------------------------------------- ##
## -------------------------------------------------------------------------- ## 
   
 
 # use all the values (will be 3000 long)
 # run four times, one for each season
   
 # spread of the values (rename these so you don't overwrite the previous vals!)
   alpha.0.sp           <- navajo_sims$alpha.0
   alpha.male.sp        <- navajo_sims$alpha.male
   
   
 # specify each season for spread - winter = 0
   alpha.season.winter.sp <- 0
   alpha.season.spring.sp <- navajo_sims$alpha.season[ , 1]
   alpha.season.summer.sp <- navajo_sims$alpha.season[ , 2]
   alpha.season.autumn.sp <- navajo_sims$alpha.season[ , 3]
   
   
 # coefficient for web   
   alpha.web.use.sp     <- c(0, navajo_sims$alpha.web[ , 1])
   
   
 # coefficients for environmental covariates (no change)
   # alpha.ndvi.sp        <- navajo_sims$alpha.ndvi
   # alpha.tmin.sp        <- navajo_sims$alpha.tmin
   # alpha.tmax.sp        <- navajo_sims$alpha.tmax
   alpha.swe.sp         <- navajo_sims$alpha.swe
   alpha.swe.winter.sp  <- navajo_sims$alpha.swe.winter
   
   
 # seasonal interactions (specify by season)
   alpha.ndvi.winter.sp <- 0
   alpha.ndvi.spring.sp <- navajo_sims$alpha.ndvi.season[ , 1]
   alpha.ndvi.summer.sp <- navajo_sims$alpha.ndvi.season[ , 2]
   alpha.ndvi.autumn.sp <- navajo_sims$alpha.ndvi.season[ , 3]
   
   
   alpha.tmin.winter.sp <- 0
   alpha.tmin.spring.sp <- navajo_sims$alpha.tmin.season[ , 1]
   alpha.tmin.summer.sp <- navajo_sims$alpha.tmin.season[ , 2]
   alpha.tmin.autumn.sp <- navajo_sims$alpha.tmin.season[ , 3]
   
   
   alpha.tmax.winter.sp <- 0
   alpha.tmax.spring.sp <- navajo_sims$alpha.tmax.season[ , 1]
   alpha.tmax.summer.sp <- navajo_sims$alpha.tmax.season[ , 2]
   alpha.tmax.autumn.sp <- navajo_sims$alpha.tmax.season[ , 3]

   
## Winter ------------------------------------------------------------------- ##   
   
   
 # female web 1
   phi.female.web1.winter <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 0 +
         alpha.season.winter.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * winter_ndvi_web1_sp +
         # alpha.tmin.sp * winter_tmin_web1_sp + 
         # alpha.tmax.sp * winter_tmax_web1_sp +
         alpha.swe.sp * winter_swe_web1_sp + 
         alpha.swe.winter.sp * winter_swewinter_web1_sp + 
         alpha.ndvi.winter.sp * winter_ndvi_web1_sp + 
         alpha.tmin.winter.sp * winter_tmin_web1_sp + 
         alpha.tmax.winter.sp * winter_tmax_web1_sp   
   )
   
   
 # male web 1 winter
   phi.male.web1.winter <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 1 +
         alpha.season.winter.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * winter_ndvi_web1_sp + 
         # alpha.tmin.sp * winter_tmin_web1_sp + 
         # alpha.tmax.sp * winter_tmax_web1_sp +
         alpha.swe.sp * winter_swe_web1_sp + 
         alpha.swe.winter.sp * winter_swewinter_web1_sp + 
         alpha.ndvi.winter.sp * winter_ndvi_web1_sp + 
         alpha.tmin.winter.sp * winter_tmin_web1_sp + 
         alpha.tmax.winter.sp * winter_tmax_web1_sp   
   )
   
   
## Spring ------------------------------------------------------------------ ##
   
   
 # female web 1
   phi.female.web1.spring <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 0 +
         alpha.season.spring.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * spring_ndvi_web1_sp + 
         # alpha.tmin.sp * spring_tmin_web1_sp + 
         # alpha.tmax.sp * spring_tmax_web1_sp +
         alpha.swe.sp * spring_swe_web1_sp + 
         alpha.swe.winter.sp * spring_swewinter_web1_sp + 
         alpha.ndvi.spring.sp * spring_ndvi_web1_sp + 
         alpha.tmin.spring.sp * spring_tmin_web1_sp + 
         alpha.tmax.spring.sp * spring_tmax_web1_sp   
   )
   
   
 # male web 1 
   phi.male.web1.spring <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 0 +
         alpha.season.spring.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * spring_ndvi_web1_sp +
         # alpha.tmin.sp * spring_tmin_web1_sp + 
         # alpha.tmax.sp * spring_tmax_web1_sp +
         alpha.swe.sp * spring_swe_web1_sp + 
         alpha.swe.winter.sp * spring_swewinter_web1_sp + 
         alpha.ndvi.spring.sp * spring_ndvi_web1_sp + 
         alpha.tmin.spring.sp * spring_tmin_web1_sp + 
         alpha.tmax.spring.sp * spring_tmax_web1_sp   
   )
   
   
## Summer ------------------------------------------------------------------- ##
   
   
 # female web 1
   phi.female.web1.summer <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 0 +
         alpha.season.summer.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * summer_ndvi_web1_sp + 
         # alpha.tmin.sp * summer_tmin_web1_sp + 
         # alpha.tmax.sp * summer_tmax_web1_sp +
         alpha.swe.sp * summer_swe_web1_sp + 
         alpha.swe.winter.sp * summer_swewinter_web1_sp + 
         alpha.ndvi.summer.sp * summer_ndvi_web1_sp + 
         alpha.tmin.summer.sp * summer_tmin_web1_sp + 
         alpha.tmax.summer.sp * summer_tmax_web1_sp   
   )
   
   
 # male web 1
   phi.male.web1.summer <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 1 +
         alpha.season.summer.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * summer_tmin_web1_sp +
         # alpha.tmin.sp * summer_tmin_web1_sp + 
         # alpha.tmax.sp * summer_tmax_web1_sp +
         alpha.swe.sp * summer_swe_web1_sp + 
         alpha.swe.winter.sp * summer_swewinter_web1_sp + 
         alpha.ndvi.summer.sp * summer_ndvi_web1_sp + 
         alpha.tmin.summer.sp * summer_tmin_web1_sp + 
         alpha.tmax.summer.sp * summer_tmax_web1_sp   
   )
   
   
## Autumn ------------------------------------------------------------------- ##
   
   
 # female web 1
   phi.female.web1.autumn <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 0 +
         alpha.season.autumn.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] + 
         # alpha.ndvi.sp * autumn_ndvi_web1_sp + 
         # alpha.tmin.sp * autumn_tmin_web1_sp + 
         # alpha.tmax.sp * autumn_tmax_web1_sp +
         alpha.swe.sp * autumn_swe_web1_sp + 
         alpha.swe.winter.sp * autumn_swewinter_web1_sp + 
         alpha.ndvi.autumn.sp * autumn_ndvi_web1_sp + 
         alpha.tmin.autumn.sp * autumn_tmin_web1_sp + 
         alpha.tmax.autumn.sp * autumn_tmax_web1_sp   
   )
   
   
 # male web 1
   phi.male.web1.autumn <- rev.logit(
      alpha.0.sp +
         alpha.male.sp * 1 +
         alpha.season.autumn.sp +            # 0 = winter (term goes to zero)
         alpha.web.use.sp[1] +
         # alpha.ndvi.sp * autumn_ndvi_web1_sp + 
         # alpha.tmin.sp * autumn_tmin_web1_sp + 
         # alpha.tmax.sp * autumn_tmax_web1_sp +
         alpha.swe.sp * autumn_swe_web1_sp + 
         alpha.swe.winter.sp * autumn_swewinter_web1_sp + 
         alpha.ndvi.autumn.sp * autumn_ndvi_web1_sp + 
         alpha.tmin.autumn.sp * autumn_tmin_web1_sp + 
         alpha.tmax.autumn.sp * autumn_tmax_web1_sp   
   )
  
   
## Violin plots for spread -------------------------------------------------- ##
   
   
 # web 1 only atm
 # combine output from above models - columns or rows?
 # female and male, winter, spring, summer, fall
   

 # rename columns
   colnames(phi.female.web1.winter) <- "FW1_WI"
   colnames(phi.female.web1.spring) <- "FW1_SP"
   colnames(phi.female.web1.summer) <- "FW1_SU"
   colnames(phi.female.web1.autumn) <- "FW1_AU"
   colnames(phi.male.web1.winter)   <- "MW1_WI"
   colnames(phi.male.web1.spring)   <- "MW1_SP"
   colnames(phi.male.web1.summer)   <- "MW1_SU"
   colnames(phi.male.web1.autumn)   <- "MW1_AU"
   
   
 # combine matrices
   web1.seasons <- cbind(phi.female.web1.winter,
                         phi.male.web1.winter,
                         phi.female.web1.spring,
                         phi.male.web1.spring,
                         phi.female.web1.summer,
                         phi.male.web1.summer,
                         phi.female.web1.autumn,
                         phi.male.web1.autumn,
                         all = TRUE)
   
   
 # make dataframe
   web1.seasons <- as.data.frame(web1.seasons)
   
    
 # remove extra column and make a date id column
   web1.seasons <- web1.seasons %>%
      select(-all) %>%
      rowid_to_column("date.id")
   
   
 # melt dataframe ()
   web1seasons <- melt(web1.seasons, id.vars = "date.id")
   
   
 # set level order 
   level_order <- factor(web1seasons$variable, level = c("FW1_WI",
                                                         "MW1_WI",
                                                         "FW1_SP",
                                                         "MW1_SP",
                                                         "FW1_SU",
                                                         "MW1_SU",
                                                         "FW1_AU",
                                                         "MW1_AU"))
   
   
## Plot spread with violin plot --------------------------------------------- ##
   

 # violin plot of web 1 male and female by seasonal category
   navajo_season_spread_plot <- ggplot(web1seasons,
                                       aes(x     = level_order, 
                                           y     = value,
                                           fill  = variable)) +
      
      geom_violin() +
      
      stat_summary(fun.data = "mean_sdl",  
                   fun.args = list(mult = 1),
                   geom = "pointrange", 
                   color = "black") + 
      
      scale_fill_viridis(discrete = TRUE) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_discrete(name   = "\nGroup",
                       breaks  = c("FW1_WI",
                                   "MW1_WI",
                                   "FW1_SP",
                                   "MW1_SP",
                                   "FW1_SU",
                                   "MW1_SU",
                                   "FW1_AU",
                                   "MW1_AU"),
                       labels = str_wrap(c("Female Winter",
                                           "Male Winter",
                                           "Female Spring",
                                           "Male Spring",
                                           "Female Summer",
                                           "Male Summer",
                                           "Female Autumn",
                                           "Male Autumn"),
                                         width = 10)) +
      
      ggtitle("Spread of Survival - Navajo\n") +
      
      theme_bw() +
      
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))

   
## Group based on male at web 1 to reduce noise ----------------------------- ##
   
   
 # combine matrices
   seasons <- cbind(phi.male.web1.winter,
                    phi.male.web1.spring,
                    phi.male.web1.summer,
                    phi.male.web1.autumn,
                    all = TRUE)
   
   
 # make dataframe
   seasons <- as.data.frame(seasons)
   
   
 # remove extra column and make a date id column
   spread_seasons <- seasons %>%
      select(-all) %>%
      rowid_to_column("date.id")
   
   
 # melt dataframe ()
   spread_seasons <- melt(spread_seasons, id.vars = "date.id")
   
   
 # set level order 
   level_order <- factor(spread_seasons$variable, level = c("MW1_WI",
                                                            "MW1_SP",
                                                            "MW1_SU",
                                                            "MW1_AU"))
   
   
## Plot assuming male at web 1 ---------------------------------------------- ##
   
   
 # violin plot of web 1 male and female by seasonal category
   navajo_season_spread_plot_mw1 <- ggplot(spread_seasons,
                                           aes(x     = level_order,
                                               y     = value,
                                               fill  = variable)) +
      
      geom_violin() +
      
      stat_summary(fun.data = "mean_sdl",  
                   fun.args = list(mult = 1),
                   geom     = "pointrange", 
                   color    = "lightgrey") +
      
      geom_boxplot(width = 0.2,
                   fill  = "lightgrey")+
      
      scale_fill_viridis(discrete = TRUE) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_discrete(name   = "\nSeason",
                       breaks  = c("MW1_WI",
                                   "MW1_SP",
                                   "MW1_SU",
                                   "MW1_AU"),
                       labels = str_wrap(c("Winter",
                                           "Spring",
                                           "Summer",
                                           "Autumn"),
                                         width = 10)) +
      
      #ggtitle("Spread of Survival - Navajo\n") +
      
      theme_bw() +
      
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold",
                                            size   = 16 ),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 12),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 12),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
   
   
## Explore effects of changing environmental variables ---------------------- ##
## -------------------------------------------------------------------------- ##
 
   
 # specify variables for new equations
   # range of values for the parameter of interest
   # use mean values for all others
   
 # use males from web 1 - with caveat in methods
 # navajo covs of interest are ndvi, tmin and tmax
 # cov_web1 dataframe to start with
   # ask range and then define vector (will always be between -1 and 1)
   
 # QUESTION 
   # when we specify range, how many values do we want/need? if we don't
   # specify, then we will just have the min and max
   
   
## Set ranges for seasons and parameters ------------------------------------ ##
## -------------------------------------------------------------------------- ##
  

## Winter ------------------------------------------------------------------- ##
   
 # specify ranges of scaled and unscaled data
   
   
 # ndvi
   winter_ndvi <- cov_joined %>%
      filter(month == 1 | month == 2 | month == 3 | month == 12) %>%
      select(ndvi, ndvi_raw) %>%
      summarise(min     = min(ndvi),
                max     = max(ndvi),
                min_raw = min(ndvi_raw),
                max_raw = max(ndvi_raw)) %>%
      mutate_all(round, 3)
   
   range_winter_ndvi     <- seq(winter_ndvi$min, winter_ndvi$max, length = 20)
   range_winter_ndvi_raw <- seq(winter_ndvi$min_raw, winter_ndvi$max_raw, length = 20)
   
   
 # tmin
   winter_tmin <- cov_joined %>%
      filter(month == 1 | month == 2 | month == 3 | month == 12) %>%
      select(tmin, tmin_raw) %>%
      summarise(min     = min(tmin),
                max     = max(tmin),
                min_raw = min(tmin_raw),
                max_raw = max(tmin_raw)) %>%
      mutate_all(round, 3)
   
   range_winter_tmin     <- seq(winter_tmin$min, winter_tmin$max, length = 20)
   range_winter_tmin_raw <- seq(winter_tmin$min_raw, winter_tmin$max_raw, length = 20)
   
   
 # tmax - note, anything with NA use na.rm (lags)
   winter_tmax <- cov_joined %>%
      filter(month == 1 | month == 2 | month == 3 | month == 12) %>%
      select(tmax, tmax_raw) %>%
      summarise(min     = min(tmax, na.rm = TRUE),
                max     = max(tmax, na.rm = TRUE),
                min_raw = min(tmax_raw, na.rm = TRUE),
                max_raw = max(tmax_raw, na.rm = TRUE)) %>%
      mutate_all(round, 3)
   
   range_winter_tmax     <- seq(winter_tmax$min, winter_tmax$max, length = 20)
   range_winter_tmax_raw <- seq(winter_tmax$min_raw, winter_tmax$max_raw, length = 20)


## Spring ------------------------------------------------------------------- ##
 
   
 # ndvi
   spring_ndvi <- cov_joined %>%
      filter(month == 4 | month == 5 | month == 6) %>%
      select(ndvi, ndvi_raw) %>%
      summarise(min     = min(ndvi),
                max     = max(ndvi),
                min_raw = min(ndvi_raw),
                max_raw = max(ndvi_raw)) %>%
      mutate_all(round, 3)
   
   range_spring_ndvi     <- seq(spring_ndvi$min, spring_ndvi$max, length = 20)
   range_spring_ndvi_raw <- seq(spring_ndvi$min_raw, spring_ndvi$max_raw, length = 20) 
   
   
 # tmin
   spring_tmin <- cov_joined %>%
      filter(month == 4 | month == 5 | month == 6) %>%
      select(tmin, tmin_raw) %>%
      summarise(min     = min(tmin),
                max     = max(tmin),
                min_raw = min(tmin_raw),
                max_raw = max(tmin_raw)) %>%
      mutate_all(round, 3)
   
   range_spring_tmin     <- seq(spring_tmin$min, spring_tmin$max, length = 20)
   range_spring_tmin_raw <- seq(spring_tmin$min_raw, spring_tmin$max_raw, length = 20)
   
   
 # tmax
   spring_tmax <- cov_joined %>%
      filter(month == 4 | month == 5 | month == 6) %>%
      select(tmax, tmax_raw) %>%
      summarise(min     = min(tmax, na.rm = TRUE),
                max     = max(tmax, na.rm = TRUE),
                min_raw = min(tmax_raw, na.rm = TRUE),
                max_raw = max(tmax_raw, na.rm = TRUE)) %>%
      mutate_all(round, 3)
   
   range_spring_tmax     <- seq(spring_tmax$min, spring_tmax$max, length = 20)
   range_spring_tmax_raw <- seq(spring_tmax$min_raw, spring_tmax$max_raw, length = 20)
   
   
## Summer ------------------------------------------------------------------- ##
   
 
 # ndvi
   summer_ndvi <- cov_joined %>%
      filter(month == 7 | month == 8 | month == 9) %>%
      select(ndvi, ndvi_raw) %>%
      summarise(min     = min(ndvi),
                max     = max(ndvi),
                min_raw = min(ndvi_raw),
                max_raw = max(ndvi_raw)) %>%
      mutate_all(round, 3)
   
   range_summer_ndvi     <- seq(summer_ndvi$min, summer_ndvi$max, length = 20)
   range_summer_ndvi_raw <- seq(summer_ndvi$min_raw, summer_ndvi$max_raw, length = 20)
   
   
 # tmin
   summer_tmin <- cov_joined %>%
      filter(month == 7 | month == 8 | month == 9) %>%
      select(tmin, tmin_raw) %>%
      summarise(min     = min(tmin),
                max     = max(tmin),
                min_raw = min(tmin_raw),
                max_raw = max(tmin_raw)) %>%
      mutate_all(round, 3)
   
   range_summer_tmin     <- seq(summer_tmin$min, summer_tmin$max, length = 20)
   range_summer_tmin_raw <- seq(summer_tmin$min_raw, summer_tmin$max_raw, length = 20)
   
   
 # tmax
   summer_tmax <- cov_joined %>%
      filter(month == 7 | month == 8 | month == 9) %>%
      select(tmax, tmax_raw) %>%
      summarise(min     = min(tmax, na.rm = TRUE),
                max     = max(tmax, na.rm = TRUE),
                min_raw = min(tmax_raw, na.rm = TRUE),
                max_raw = max(tmax_raw, na.rm = TRUE)) %>%
      mutate_all(round, 3)
   
   range_summer_tmax     <- seq(summer_tmax$min, summer_tmax$max, length = 20)
   range_summer_tmax_raw <- seq(summer_tmax$min_raw, summer_tmax$max_raw, length = 20) 
    
   
## Autumn ------------------------------------------------------------------- ##

   
 # ndvi
   autumn_ndvi <- cov_joined %>%
      filter(month == 10 | month == 11) %>%
      select(ndvi, ndvi_raw) %>%
      summarise(min     = min(ndvi),
                max     = max(ndvi),
                min_raw = min(ndvi_raw),
                max_raw = max(ndvi_raw)) %>%
      mutate_all(round, 3)
   
   range_autumn_ndvi     <- seq(autumn_ndvi$min, autumn_ndvi$max, length = 20)
   range_autumn_ndvi_raw <- seq(autumn_ndvi$min_raw, autumn_ndvi$max_raw, length = 20)
   
   
 # tmin
   autumn_tmin <- cov_joined %>%
      filter(month == 10 | month == 11) %>%
      select(tmin, tmin_raw) %>%
      summarise(min     = min(tmin),
                max     = max(tmin),
                min_raw = min(tmin_raw),
                max_raw = max(tmin_raw)) %>%
      mutate_all(round, 3)
   
   range_autumn_tmin     <- seq(autumn_tmin$min, autumn_tmin$max, length = 20)
   range_autumn_tmin_raw <- seq(autumn_tmin$min_raw, autumn_tmin$max_raw, length = 20)
   
   
 # tmax
   autumn_tmax <- cov_joined %>%
      filter(month == 10 | month == 11) %>%
      select(tmax, tmax_raw) %>%
      summarise(min     = min(tmax, na.rm = TRUE),
                max     = max(tmax, na.rm = TRUE),
                min_raw = min(tmax_raw, na.rm = TRUE),
                max_raw = max(tmax_raw, na.rm = TRUE)) %>%
      mutate_all(round, 3)
   
   range_autumn_tmax     <- seq(autumn_tmax$min, autumn_tmax$max, length = 20)
   range_autumn_tmax_raw <- seq(autumn_tmax$min_raw, autumn_tmax$max_raw, length = 20)   
   
   
## Specify effect of ndvi, tmin, and tmax by season ------------------------- ##
## -------------------------------------------------------------------------- ##
   
   
## Winter ------------------------------------------------------------------- ##
   

 # assumes males from web 1, season winter = 1 here although in the model run
 # it's specified as 0
   
   # NOTE - does not work unless alpha.tmin.season.use[1] is specified BUT that
   # means the entire length of range_winter_tmin is multiplied by 0
   # but we should still have 15 values for phi given the
   # alpha.tmin * range_winter_tmin line of the equation
   # also, if we have to specify alpha.tmin.season.use here, why not above?
   
   # maybe not working because other params are length of 4 - set season for all
   # REMEMBER - defined above, first item in vector [1] = 0 for winter
   
   
 # effect of ndvi
   phi.eff.ndvi.winter <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[1] + 
         alpha.web.use[1] + 
         # alpha.ndvi * range_winter_ndvi + 
         # alpha.tmin * mean.tmin.w1[1] + 
         # alpha.tmax * mean.tmax.w1[1] +
         alpha.swe * mean.swe.w1[1] + 
         alpha.swe.winter * mean.swe.winter.w1[1] + 
         alpha.ndvi.season.use[1] * range_winter_ndvi +
         alpha.tmin.season.use[1] * mean.tmin.w1[1] +
         alpha.tmax.season.use[1] * mean.tmax.w1[1]
      
   )
   
 # effect of tmin
   phi.eff.tmin.winter <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[1] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[1] +
         # alpha.tmin * range_winter_tmin + 
         # alpha.tmax * mean.tmax.w1[1] +
         alpha.swe * mean.swe.w1[1] + 
         alpha.swe.winter * mean.swe.winter.w1[1] + 
         alpha.ndvi.season.use[1] * mean.ndvi.w1[1] +
         alpha.tmin.season.use[1] * range_winter_tmin +
         alpha.tmax.season.use[1] * mean.tmax.w1[1]
      
      )
   
   
 # effect of tmax
   phi.eff.tmax.winter <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[1] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[1] + 
         # alpha.tmin * mean.tmin.w1[1] + 
         # alpha.tmax * range_winter_tmax +
         alpha.swe * mean.swe.w1[1] + 
         alpha.swe.winter * mean.swe.winter.w1[1] + 
         alpha.ndvi.season.use[1] * mean.ndvi.w1[1] +
         alpha.tmin.season.use[1] * mean.tmin.w1[1] +
         alpha.tmax.season.use[1] * range_winter_tmax
      
   )
   
## Spring ------------------------------------------------------------------- ##
   
   
 # effect of ndvi
   phi.eff.ndvi.spring <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[2] + 
         alpha.web.use[1] + 
         # alpha.ndvi * range_spring_ndvi + 
         # alpha.tmin * mean.tmin.w1[2] +
         # alpha.tmax * mean.tmax.w1[2] +
         alpha.swe * mean.swe.w1[2] + 
         alpha.swe.winter * mean.swe.winter.w1[2] + 
         alpha.ndvi.season.use[2] * range_spring_ndvi +
         alpha.tmin.season.use[2] * mean.tmin.w1[2] +
         alpha.tmax.season.use[2] * mean.tmax.w1[2]
   )
   
   
 # effect of tmin
   phi.eff.tmin.spring <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[2] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[2] +
         # alpha.tmin * range_spring_tmin + 
         # alpha.tmax * mean.tmax.w1[2] +
         alpha.swe * mean.swe.w1[2] + 
         alpha.swe.winter * mean.swe.winter.w1[2] + 
         alpha.ndvi.season.use[2] * mean.ndvi.w1[2] +
         alpha.tmin.season.use[2] * range_spring_tmin +
         alpha.tmax.season.use[2] * mean.tmax.w1[2]
      
   )
   
   
   # effect of tmax
   phi.eff.tmax.spring <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[2] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[2] + 
         # alpha.tmin * mean.tmin.w1[2] + 
         # alpha.tmax * range_spring_tmax +
         alpha.swe * mean.swe.w1[2] + 
         alpha.swe.winter * mean.swe.winter.w1[2] +
         alpha.ndvi.season.use[2] * mean.ndvi.w1[2] +
         alpha.tmin.season.use[2] * mean.tmin.w1[2] +
         alpha.tmax.season.use[2] * range_spring_tmax
   )
   
   
## Summer ------------------------------------------------------------------- ##
   
   
 # effect of ndvi
   phi.eff.ndvi.summer <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[3] + 
         alpha.web.use[1] + 
         # alpha.ndvi * range_summer_ndvi + 
         # alpha.tmin * mean.tmin.w1[3] +
         # alpha.tmax * mean.tmax.w1[3] +
         alpha.swe * mean.swe.w1[3] + 
         alpha.swe.winter * mean.swe.winter.w1[1] + 
         alpha.ndvi.season.use[3] * range_summer_ndvi +
         alpha.tmin.season.use[3] * mean.tmin.w1[3] +
         alpha.tmax.season.use[3] * mean.tmax.w1[3]
   )
   
   
 # effect of tmin
   phi.eff.tmin.summer <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[3] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[3] +
         # alpha.tmin * range_summer_tmin + 
         # alpha.tmax * mean.tmax.w1[3] +
         alpha.swe * mean.swe.w1[3] + 
         alpha.swe.winter * mean.swe.winter.w1[1] + 
         alpha.ndvi.season.use[3] * mean.ndvi.w1[3] + 
         alpha.tmin.season.use[3] * range_summer_tmin +
         alpha.tmax.season.use[3] * mean.tmax.w1[3]
   )
   
   
 # effect of tmax
   phi.eff.tmax.summer <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[3] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[3] + 
         # alpha.tmin * mean.tmin.w1[3] + 
         # alpha.tmax * range_summer_tmax +
         alpha.swe * mean.swe.w1[3] + 
         alpha.swe.winter * mean.swe.winter.w1[3] + 
         alpha.ndvi.season.use[3] * mean.ndvi.w1[3] +
         alpha.tmin.season.use[3] * mean.tmin.w1[3] +
         alpha.tmax.season.use[3] * range_summer_tmax
   )
   
   
## Autumn ------------------------------------------------------------------- ##

   
 # effect of ndvi
   phi.eff.ndvi.autumn <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[4] + 
         alpha.web.use[1] + 
         # alpha.ndvi * range_autumn_ndvi + 
         # alpha.tmin * mean.tmin.w1[4] +
         # alpha.tmax * mean.tmax.w1[4] +
         alpha.swe * mean.swe.w1[4] + 
         alpha.swe.winter * mean.swe.winter.w1[4] + 
         alpha.ndvi.season.use[4] * range_autumn_ndvi +
         alpha.tmin.season.use[4] * mean.tmin.w1[4] +
         alpha.tmax.season.use[4] * mean.tmax.w1[4]
   )
   
   
 # effect of tmin
   phi.eff.tmin.autumn <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[4] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[4] +
         # alpha.tmin * range_autumn_tmin + 
         # alpha.tmax * mean.tmax.w1[4] +
         alpha.swe * mean.swe.w1[4] + 
         alpha.swe.winter * mean.swe.winter.w1[4] + 
         alpha.ndvi.season.use[4] * mean.ndvi.w1[4] +
         alpha.tmin.season.use[4] * range_autumn_tmin +
         alpha.tmax.season.use[4] * mean.tmax.w1[4]
   )
   
   
 # effect of tmax
   phi.eff.tmax.autumn <- rev.logit(
      alpha.0 +
         alpha.male * 1 +                      
         alpha.season.use[4] + 
         alpha.web.use[1] + 
         # alpha.ndvi * mean.ndvi.w1[4] + 
         # alpha.tmin * mean.tmin.w1[4] + 
         # alpha.tmax * range_autumn_tmax +
         alpha.swe * mean.swe.w1[4] + 
         alpha.swe.winter * mean.swe.winter.w1[4] + 
         alpha.ndvi.season.use[4] * mean.ndvi.w1[4] +
         alpha.tmin.season.use[4] * mean.tmin.w1[4] +
         alpha.tmax.season.use[4]* range_autumn_tmax
   )
   
## Plot effect of tmin and tmax on phi -------------------------------------- ##
## -------------------------------------------------------------------------- ##
  
   
## UPDATE - facet the sig. stuff -------------------------------------------- ##
   

 # ndvi
   phi.eff.ndvi <- ggplot() +
      geom_line(aes(x = range_autumn_ndvi_raw , y = phi.eff.ndvi.autumn,
                    color = "Autumn"), size = 2.0) +
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      scale_x_continuous(name    = "\nNormalized difference vegetation index \n ndvi0",
                         breaks  = seq(0, 0.4, 0.1),
                         limits  = c(0, 0.4)) +
      scale_colour_manual(breaks = "Autumn",
                          values = "#FDE725FF") + 
      theme(legend.title = element_blank(),
            legend.position="bottom",
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            axis.title       = element_text(face   = "bold", size   = 12),
            axis.text.y      = element_text(colour = "black", size   = 12),
            axis.line        = element_line(color  = "black", size   = 0.5))
   

   
## Plot phi with unscaled covariates ---------------------------------------- ##
   

 # ndvi
   phi.eff.ndvi.plot <- ggplot() +
      
      # create lines - this is a total hack :(
      geom_line(aes(x = range_winter_ndvi_raw , y = phi.eff.ndvi.winter,
                    color = "Winter"), size = 1.5) +
      geom_line(aes(x = range_spring_ndvi_raw , y = phi.eff.ndvi.spring,
                    color = "Spring"), size = 1.5) +
      geom_line(aes(x = range_summer_ndvi_raw , y = phi.eff.ndvi.summer,
                    color = "Summer"), size = 1.5) +
      geom_line(aes(x = range_autumn_ndvi_raw , y = phi.eff.ndvi.autumn,
                    color = "Autumn"), size = 1.5) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_continuous(name    = "\nNDVI",
                         breaks  = seq(-1, 1, 0.2),
                         limits  = c(-1, 1)) +
      
      labs(title   = "Navajo",
           subtitle = "Average monthly NDVI") +
      
      scale_colour_manual("",
                          breaks = c("Winter",
                                     "Spring",
                                     "Summer",
                                     "Autumn"),
                          values = c("#440154FF",
                                     "#31688EFF",
                                     "#35B779FF",
                                     "#FDE725FF")) + 
      
      theme_bw() +
      
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            plot.subtitle    = element_text(color  = "gray28"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
   
   
   
 # tmin
   phi.eff.tmin.plot <- ggplot() +
      
      # create lines - this is a total hack :(
      geom_line(aes(x = range_winter_tmin_raw , y = phi.eff.tmin.winter,
                    color = "Winter"), size = 1.5) +
      geom_line(aes(x = range_spring_tmin_raw , y = phi.eff.tmin.spring,
                    color = "Spring"), size = 1.5) +
      geom_line(aes(x = range_summer_tmin_raw , y = phi.eff.tmin.summer,
                    color = "Summer"), size = 1.5) +
      geom_line(aes(x = range_autumn_tmin_raw , y = phi.eff.tmin.autumn,
                    color = "Autumn"), size = 1.5) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_continuous(name    = "\nTemperature (\u00B0C)",
                         breaks  = seq(-14, 14, 2),
                         limits  = c(-14, 14)) +

      labs(title    = "Navajo",
           subtitle = "Average minimum temperature (\u00B0C)") +
      
      scale_colour_manual("",
                          breaks = c("Winter",
                                     "Spring",
                                     "Summer",
                                     "Autumn"),
                          values = c("#440154FF",
                                     "#31688EFF",
                                     "#35B779FF",
                                     "#FDE725FF")) + 
      
      theme_bw() +
      
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            plot.subtitle    = element_text(color  = "gray28"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
   

 # tmax
   phi.eff.tmax.plot <- ggplot() +
      
      # create lines - this is a total hack :(
      geom_line(aes(x = range_winter_tmax_raw , y = phi.eff.tmax.winter,
                    color = "Winter"), size = 1.5) +
      geom_line(aes(x = range_spring_tmax_raw , y = phi.eff.tmax.spring,
                    color = "Spring"), size = 1.5) +
      geom_line(aes(x = range_summer_tmax_raw , y = phi.eff.tmax.summer,
                    color = "Summer"), size = 1.5) +
      geom_line(aes(x = range_autumn_tmax_raw , y = phi.eff.tmax.autumn,
                    color = "Autumn"), size = 1.5) +
      
      scale_y_continuous(name   = "Survival\n",
                         breaks = seq(0, 1.0, 0.1),
                         limits = c(0, 1.0)) +
      
      scale_x_continuous(name    = "\nTemperature (\u00B0C)",
                         breaks  = seq(-0, 35, 5),
                         limits  = c(0, 35)) +
      
      labs(title    = "Navajo",
           subtitle = "Average maximum temperature (\u00B0C)") +
      
      scale_colour_manual("",
                          breaks = c("Winter",
                                     "Spring",
                                     "Summer",
                                     "Autumn"),
                          values = c("#440154FF",
                                     "#31688EFF",
                                     "#35B779FF",
                                     "#FDE725FF")) + 
      
      theme_bw() +
      
      theme(legend.title = element_blank(),
            panel.grid.major = element_line(color = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border     = element_blank(), 
            panel.background = element_blank(),
            plot.title       = element_text(size   = 14, 
                                            family = "Tahoma", 
                                            face   = "bold"),
            plot.subtitle    = element_text(color  = "gray28"),
            text             = element_text(family = "Tahoma"),
            axis.title       = element_text(face   = "bold"),
            axis.text.x      = element_text(color  = "black", 
                                            size   = 11),
            axis.text.y      = element_text(colour = "black", 
                                            size   = 9),
            axis.line        = element_line(size   = 0.5, 
                                            color  = "black"))
  
   
## Look at distribution of covariates in relation to phi -------------------- ##
## -------------------------------------------------------------------------- ##
   
   
 # make object with vector for seasons
   w <- which(covariate.data$temporal.covariates$season == 1)
   s <- which(covariate.data$temporal.covariates$season == 3)
   a <- which(covariate.data$temporal.covariates$season == 4)
   
   
 # now subset phi by season - use male web 1 index 149
   phi.w <- navajo_model$BUGSoutput$sims.list$phi.male.web1[149, w]
   phi.s <- navajo_model$BUGSoutput$sims.list$phi.male.web1[149, s]
   phi.a <- navajo_model$BUGSoutput$sims.list$phi.male.web1[149, a]
   
   
 # subset ndvi by season
   winter <- covariate.data$temporal.covariates$ndvi_0.webnavajo.1[
      covariate.data$temporal.covariates$season %in% "1"]
   
   summer <- covariate.data$temporal.covariates$ndvi_0.webnavajo.1[
      covariate.data$temporal.covariates$season %in% "3"]
   
   autumn <- unscale_covariates %>%
      filter(site == "navajo") %>%
      filter(web == 1) %>%
      filter(month == 10 | month == 11) %>%
      select(ndvi)
      
      #covariate.data$temporal.covariates$ndvi_0.webnavajo.1[
      #covariate.data$temporal.covariates$season %in% "4"]
   
 # can't plot like this since they are different lengths
   plot(winter, phi.w, xlim = c(-6, 6), ylim = c(0, 1))
   plot(summer, phi.s, xlim = c(-6, 6), ylim = c(0, 1))
   plot(autumn, phi.a, xlim = c(-2, 2), ylim = c(0, 1))
   
   
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##