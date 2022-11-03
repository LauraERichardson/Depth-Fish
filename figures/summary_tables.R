library(tidyverse)
library(brms)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############ 

# Summary tables S4, S5, S6, S11
sink("summary_tables.txt", append = F)
summaries <- lapply(1:length(models), function(i){
  model <- models[i]
  cat('\n\nModel: ', model, ' - TABLE  S4 summary\n\n')
  print(summary(get(model)))
  cat('\n\nModel: ', model, ' - TABLE  S5 directional P\n\n')
  print(bayestestR::p_direction(get(model), effects="fixed"))
  cat('\n\nModel: ', model, ' - TABLE  S6 Bayes R2\n\n')
  performance::r2_bayes(get(model))
  cat('\n\nModel: ', model, ' - TABLE  S11 Variance partitioning\n\n')
  sds <- posterior_samples(get(model), pars = 'sd')
  print(length(which(sds$sd_ECOREGION__Intercept > sds$sd_ISLAND__Intercept))/length(sds$sd_ECOREGION__Intercept))
  print(length(which(sds$sd_ECOREGION__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ECOREGION__Intercept))
  print(length(which(sds$sd_ISLAND__Intercept > sds$sd_SITE__Intercept))/length(sds$sd_ISLAND__Intercept))
  if(hurd[i]){
    cat('\n\nModel: ', model, ' - TABLE  S11 Variance partitioning - Hurdle component\n\n')
    length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_ISLAND__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
    length(which(sds$sd_ECOREGION__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ECOREGION__hu_Intercept)
    length(which(sds$sd_ISLAND__hu_Intercept > sds$sd_SITE__hu_Intercept))/length(sds$sd_ISLAND__hu_Intercept)
    
  }
})

###################################################################

# End
