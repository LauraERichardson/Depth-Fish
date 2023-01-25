library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############## 

# Figure X: depth effect by island: ggplot version

ppd <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  newdat <- data.frame(distinct(set %>% select(ECOREGION, ISLAND, POP_STATUS)))
  newdat$trophic_group <- factor(nms[i], levels = nms)
  
  re <- ranef(get(models[i]))
  
  re$ISLAND[,,'DEPTH_c'] %>% as.data.frame() %>% mutate(ISLAND=rownames(.)) %>% inner_join(newdat)
  
}, mc.cores=5) %>% bind_rows()


ppd <- ppd %>% group_by(trophic_group, ECOREGION) %>% arrange(ECOREGION, POP_STATUS, Estimate) 
ppd$ISLAND <- forcats::fct_inorder(ppd$ISLAND)
ppds <- ppd %>%
  group_by(trophic_group) %>%
  mutate(ymins = min(Q2.5)-0.049,
         ymaxs = max(Q97.5)+0.049) %>%
  group_by(trophic_group,ECOREGION) %>%
  summarise(xmins = min(as.numeric(ISLAND))-0.499,
            xmaxs = max(as.numeric(ISLAND))+0.499,
            ymins = unique(ymins),
            ymaxs = unique(ymaxs))

ggplot() + 
  geom_point(aes(x=ISLAND, y=Estimate, col=trophic_group), data=ppd) +
  geom_linerange(aes(x=ISLAND, y=Estimate, ymax=Q97.5, ymin=Q2.5, col=trophic_group,size=POP_STATUS), data=ppd) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_colour_manual('',values = mycols, guide='none') +
  geom_rect(aes(xmin=xmins, xmax=xmaxs, ymin=ymins, ymax=ymaxs, fill=trophic_group, col= trophic_group,alpha=ECOREGION), data=ppds) +
  scale_alpha_discrete('EcoRegion',range=c(0.05,0.6),breaks=function(x) rev(x)) +
  #scale_linewidth_discrete('',range=c(0.5,1.1), guide='none') +
  scale_size_discrete('',range=c(0.6,1.1), guide='none') +
  facet_wrap(~trophic_group, scales='free', nrow = 1) +
  geom_hline(yintercept = 0, linetype=2) +
  coord_flip() +
  xlab('Island') + 
  ylab("Depth effect") +
  cowplot::theme_cowplot() + 
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text.x = element_text(face = 'bold')
  )

ggsave('Figure_Isl_ER_depthfx.pdf', device="pdf", width = 20, height = 12, units = 'in',dpi = 300, bg = "white")
#ggsave('Figure_Isl_ER_depthfx.png', device="png", width = 20, height = 12, units = 'in',dpi = 300, bg = "white")



###################################################################

# End
