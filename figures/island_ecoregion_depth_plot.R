library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############## 

# Figure X: depth effect by island: ggplot version

ppd_l <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  newdat <- data.frame(distinct(set %>% select(ECOREGION, ISLAND, POP_STATUS)) %>% 
                         mutate(DEPTH_c=1,
                                SITE_SLOPE_400m_c=0))
  newdat$trophic_group <- factor(nms[i], levels = nms)
  
  nd1 <- newdat %>% 
    add_epred_draws(get(models[i]), 
                    re_formula = ~(1+DEPTH_c||ISLAND), seed = 123, dpar=TRUE) %>%
    mutate(IFX = 1)
  
  nd2 <- newdat %>% 
    add_epred_draws(get(models[i]), 
                    re_formula = NA, seed = 123, dpar=TRUE)%>%
    mutate(IFX = 0)
  
  bind_rows(nd1, nd2) %>% group_by(trophic_group, ECOREGION, ISLAND, POP_STATUS,.row, .draw) %>% summarise(.epred = .epred[1]/.epred[2])
  
}, mc.cores=5) %>% bind_rows()


ppd <- ppd_l %>% median_qi(.epred) %>% group_by(trophic_group, ECOREGION) %>% arrange(ECOREGION, POP_STATUS, .epred) 
ppd$ISLAND <- forcats::fct_inorder(ppd$ISLAND)
ppds <- ppd %>%
  group_by(trophic_group) %>%
  mutate(ymins = min(.lower)-0.049,
         ymaxs = max(.upper)+0.049) %>%
  group_by(trophic_group,ECOREGION) %>%
  summarise(xmins = min(as.numeric(ISLAND))-0.499,
            xmaxs = max(as.numeric(ISLAND))+0.499,
            ymins = unique(ymins),
            ymaxs = unique(ymaxs))

ggplot() + 
  geom_point(aes(x=ISLAND, y=.epred, col=trophic_group), data=ppd) +
  geom_linerange(aes(x=ISLAND, y=.epred, ymax=.upper, ymin=.lower, col=trophic_group,size=POP_STATUS), data=ppd) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_colour_manual('',values = mycols, guide='none') +
  geom_rect(aes(xmin=xmins, xmax=xmaxs, ymin=ymins, ymax=ymaxs, fill=trophic_group, col= trophic_group,alpha=ECOREGION), data=ppds) +
  scale_alpha_discrete('EcoRegion',range=c(0.05,0.6),breaks=function(x) rev(x)) +
  #scale_linewidth_discrete('',range=c(0.5,1.1), guide='none') +
  scale_size_discrete('',range=c(0.6,1.1), guide='none') +
  facet_wrap(~trophic_group, scales='free', nrow = 1) +
  geom_hline(yintercept = 1, linetype=2) +
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
