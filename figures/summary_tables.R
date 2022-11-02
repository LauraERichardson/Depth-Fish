library(brms)
library(bayestestR)
library(performance)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############ 

# Summary tables S4, S5, S6, S11

summaries <- lapply(models, function(models){
  list(tableS4 = summary(models),
       tableS5 = p_direction(models, effects="fixed"),
       tableS6 = performance::r2_bayes(models))#,
       #tableS11 = sds <- posterior_samples(models, pars = 'sd'),
       #length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept),
       #length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept),
       #length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept)) 
       }
  )
  
write_csv(as.data.frame(summaries), 'summary_tables.csv')

###################################################################

# End







# Table S11: Probabilities of differing variation in fish biomass across hierarchical spatial scales: ecoregion, island, and site

# Total fish biomass

sds <- posterior_samples(fish_mod, pars = 'sd')
#Probability that ecoregion is different to island, and site:
length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept) 
length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept) 
#Probability that island is different to site:
length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept) 

# Primary consumer biomass

sds <- posterior_samples(PRIM_mod, pars = 'sd')
#Probability that ecoregion is different to island, and site:
length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept)
length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept)
#Probability that island is different to site:
length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept)

# Planktivore biomass

sds <- posterior_samples(PLANK_mod, pars = 'sd')
#Probability that ecoregion is different to island, and site:
length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept)
length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept) 
length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_ISLAND__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
#Probability that island is different to site:
length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept)
length(which(sds$sd_ISLAND__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ISLAND__hu_Intercept)

# Secondary consumer biomass

sds <- posterior_samples(SEC_mod, pars = 'sd')
#Probability that ecoregion is different to island, and site:
length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept)
length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept)
#Probability that island is different to site:
length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept)

# Piscivore biomass

sds <- posterior_samples(PISC_mod, pars = 'sd')
#Probability that ecoregion is different to island, and site:
length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept)
length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept)
length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_ISLAND__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
#Probability that island is different to site:
length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept)
length(which(sds$sd_ISLAND__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ISLAND__hu_Intercept)
