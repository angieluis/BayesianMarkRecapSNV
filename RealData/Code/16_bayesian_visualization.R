## Bayesian visualization with various packages ----------------------------- ##


 # packages
   library(coda)
   library(hexbin)
   library(bayesplot)
   library(MCMCvis)
   library(tidybayes)
   

zuni_fit <- as.mcmc(zuni_model)


# read model output into coda format
navajo_fit <- as.mcmc(navajo_model)


mcmc_areas(navajo_fit, pars = c("alpha.0"),prob = 0.90) 

mcmc_areas(zuni_fit, pars = c("alpha.tmin.season[1]"),prob = 0.90) 

mcmc_scatter(navajo_fit, pars = c("alpha.ndvi.season[3]", "alpha.0"),
             size = 1.5, alpha = 0.5)

if (requireNamespace("hexbin", quietly = TRUE)) {
  mcmc_hex(navajo_fit, pars = c("alpha.ndvi.season[3]", "alpha.0"))
}


color_scheme_set("darkgray")
if (requireNamespace("hexbin", quietly = TRUE)) {
  mcmc_hex(navajo_fit, pars = c("alpha.season[1]", "alpha.0"))
}


exp(x)/(1+exp(x))


