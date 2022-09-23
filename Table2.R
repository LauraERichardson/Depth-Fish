library(tidyverse)
library(brms)
library(tidybayes)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############ # Figure 3: Probability of proportional increase  (across 0-30m depth; at unpop islands; with slope held constant)

# set up data frame to predict from the posterior across depths 0 to 30m in 10m jumps, with slope held constant (at mean)
# calculate the changes 

pp <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  newdat <- data.frame(expand.grid("POP_STATUS"=levels(set$POP_STATUS), 
                                   "DEPTH"=seq(0,30, by=10), 
                                   "OBS_YEAR"=levels(set$OBS_YEAR)[1], 
                                   "SITE_SLOPE_400m_c"=mean(set$SITE_SLOPE_400m_c,na.rm=TRUE), 
                                   "ISLAND"='FOO', 
                                   "ECOREGION"='FOO', 
                                   "SITE"='FOO',
                                   "DIVER"='FOO'))
  newdat$DEPTH_c <- (newdat$DEPTH - mean(set$DEPTH))/sd(set$DEPTH)
  newdat$trophic_group <- factor(nms[i], levels = nms)
  
  newdat %>% 
    add_epred_draws(get(models[i]), 
                    re_formula = NA, seed = 123, dpar=TRUE)
  
}, mc.cores = 5) %>% bind_rows()


pp %>% group_by(trophic_group, .draw, POP_STATUS) %>%
  select(-OBS_YEAR,-SITE_SLOPE_400m_c,-ISLAND,-ECOREGION,-SITE,- DIVER,-.chain,-.iteration) %>%
  arrange(trophic_group, POP_STATUS, .draw, DEPTH) %>%
  mutate(
    hu = ifelse(is.na(hu),0,hu),
    expect = (1-hu)*mu,
    rat = expect/lag(expect)) %>%
  filter(!is.na(rat)) %>%
  arrange(trophic_group, .draw, DEPTH,POP_STATUS) %>%
  group_by(trophic_group, .draw, DEPTH) %>%
  summarise(rat_pop = rat[2]/rat[1]) %>% # populated over unpopulated
  group_by(trophic_group, DEPTH) %>%
  summarise(prob=mean(rat_pop>1)) %>%  # probability that fish biomass at populated islands increases faster than at unpopulated islands
  pivot_wider(names_from = trophic_group, values_from = prob)
  

