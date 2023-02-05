library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

##########  

pphus <- parallel::mclapply(1:length(models), function(i) {
 
  post_int <- as_draws_df(get(models[i]),'b_(.*_)?Intercept',regex = T) %>% 
    mutate(Intercept = exp(b_Intercept)) %>% 
    mutate(across(contains('b_hu_Intercept'),~(1-inv_logit(.x))*Intercept))
  
  iter <- nrow(post_int)
   
  pos_prior <- get(models[i])$prior[get(models[i])$prior$class=='Intercept','prior'][1]
  pm <- as.numeric(gsub(".*\\(([0-9].[0-9]*),.*","\\1",pos_prior))
  psd <- as.numeric(gsub(".*\\([0-9].[0-9]*,([0-9]).*","\\1",pos_prior))
  
  if(input_frame$HURDLE[i]){
    post_int$priors <- inv_logit(rlogis(iter, -2, 0.5)) * exp(rnorm(iter, pm, psd))
  } else {
    post_int$priors <- exp(rnorm(iter, pm, psd))
  }
  
  post_int$trophic_group = factor(nms[i], levels = nms)
  post_int
  
}, mc.cores = 5) %>% bind_rows()

ggplot(pphus) +
  geom_histogram(aes(x=ifelse(!is.na(b_hu_Intercept),b_hu_Intercept,Intercept), y = after_stat(density)), bins = 30, fill='orange2') + 
  geom_density(aes(x=priors, y = after_stat(density)), col='skyblue') + 
  scale_x_log10() +
  facet_wrap(~trophic_group, scales='free', ncol=3) + 
  ylab('Density') +
  xlab('Biomass density (kg/ha)') +
  cowplot::theme_cowplot() 


ggsave('Prior_post_dens.pdf',device="pdf", width = 10, height = 6, units = 'in',dpi = 300)

