library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

if(!require(ggridges)) install.packages('ggridges')

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############ # Figure 3: Probability of proportional increase  (across 0-30m depth; at unpop islands; with slope held constant)


pp <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  newdat <- data.frame(expand.grid("POP_STATUS"=levels(set$POP_STATUS), 
                                   "DEPTH"=seq(0,30, by=10), 
                                   "SITE_SLOPE_400m_c"=mean(set$SITE_SLOPE_400m_c,na.rm=TRUE),
                                   "ISLAND"='FOO', 
                                   "ECOREGION"='FOO', 
                                   "SITE"='FOO',
                                   "DIVER"='FOO',
                                   "OBS_YEAR"='FOO'))
  
  newdat$DEPTH_c <- (newdat$DEPTH - mean(set$DEPTH))/sd(set$DEPTH)
  newdat$trophic_group <- factor(nms[i], levels = nms)
  
  newdat %>% 
    add_epred_draws(get(models[i]), 
                    re_formula = NA, seed = 123, dpar=TRUE)
  
}, mc.cores = 5) %>% bind_rows()


pp_bd <- pp %>% group_by(trophic_group, .draw, POP_STATUS) %>%
  select(-OBS_YEAR,-SITE_SLOPE_400m_c,-ISLAND,-ECOREGION,-SITE,- DIVER,-.chain,-.iteration) %>%
  arrange(trophic_group, POP_STATUS, .draw, DEPTH) %>%
  mutate(
    hu = ifelse(is.na(hu),0,hu),
    expect = (1-hu)*mu,
    rat = expect - lag(expect)) %>%
  filter(!is.na(rat)) %>%
  arrange(trophic_group, .draw, DEPTH,POP_STATUS) %>%
  group_by(trophic_group, .draw, DEPTH) %>%
  summarise(rat_pop = rat[1] - rat[2])

g1 <- pp_bd %>% 
  ggplot() + 
  geom_histogram(aes(x=rat_pop, fill=trophic_group, after_stat(density), group=factor(DEPTH), alpha=as.factor(DEPTH)),bins = 30) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_alpha_discrete('',labels = c('0m-10m','10m-20m','20m-30m')) +
  facet_wrap(~trophic_group, scales='free', ncol = 1) +
  xlab('Zonation ratio') + 
  ylab('') +
  geom_vline(xintercept=1, linetype=2) +
  cowplot::theme_cowplot() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text()
  )

pp_bdp <- pp %>% group_by(trophic_group, .draw, POP_STATUS) %>%
  select(-OBS_YEAR,-SITE_SLOPE_400m_c,-ISLAND,-ECOREGION,-SITE,- DIVER,-.chain,-.iteration) %>%
  arrange(trophic_group, POP_STATUS, .draw, DEPTH) %>%
  mutate(
    hu = ifelse(is.na(hu),0,hu),
    expect = (1-hu)*mu,
    rat = expect - lag(expect)) %>%
  filter(!is.na(rat)) %>%
  group_by(trophic_group) %>% 
  filter(rat <= quantile(rat,0.99))%>% 
  mutate(rat = rat - 0)

g2 <- ggplot() + 
  geom_density_ridges2(
    aes(x=rat*100, 
        y = POP_STATUS,
        fill=trophic_group,
        col=trophic_group
    ),alpha=0.4,data=pp_bdp) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_colour_manual('',values = mycols, guide='none') +
  scale_alpha_discrete('') +
  facet_wrap(trophic_group~DEPTH, scales='free', ncol=3) + 
  ylab('Density') +
  xlab('Absolute change with depth bin') +
  geom_vline(xintercept=0, linetype=2) +
  cowplot::theme_cowplot() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(colour = "black", fill = "white"),
    strip.background = element_rect(colour = "black", fill = "white")
  )

g2+g1 + patchwork::plot_layout(widths=c(0.75,0.25))

ggsave('Figure3_alt.png',width = 7, height = 6, units = 'in',dpi = 150)

pp_bdp2 <- pp %>% group_by(trophic_group, .draw, POP_STATUS) %>%
  select(-OBS_YEAR,-SITE_SLOPE_400m_c,-ISLAND,-ECOREGION,-SITE,- DIVER,-.chain,-.iteration) %>%
  arrange(trophic_group, POP_STATUS, .draw, DEPTH) %>%
  mutate(
    hu = ifelse(is.na(hu),0,hu),
    expect = (1-hu)*mu,
    rat = expect/lag(expect)) %>%
  filter(!is.na(rat)) %>%
  group_by(trophic_group, DEPTH, POP_STATUS) %>%
  summarise('10%' = mean(rat>1.1),
            '20%' = mean(rat>1.2),
            '50%' = mean(rat>1.5),
            '100%'   = mean(rat>2),
            '200%'   = mean(rat>3)) %>%
  pivot_longer(cols = contains('%'), names_to = 'increase', values_to = 'ratio')

pp_bdp2$increase <- factor(pp_bdp2$increase, levels = c('10%', '20%', '50%', '100%', '200%'))

pp_bdp2 %>% 
  filter(DEPTH>0) %>%
  ggplot() + 
  geom_point(aes(x=DEPTH, y=ratio, col=trophic_group)) +
  geom_line(aes(x=DEPTH, y=ratio, col=trophic_group)) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_colour_manual('',values = mycols, guide='none') +
  scale_alpha_discrete('') +
  facet_grid(POP_STATUS~increase) +
  xlab('Depth (m)') + 
  ylab('Probability') +
  geom_hline(yintercept=0.75, linetype=2) +
  geom_hline(yintercept=0.95, linetype=3) +
  geom_hline(yintercept=0) +
  cowplot::theme_cowplot() + 
  coord_cartesian(ylim=c(00.05,1),expand = 0.05,clip = 'off') +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )

ggsave('Figure3_gg.png',width = 7, height = 6, units = 'in',dpi = 150)

png('Figure3.png',width = 7, height = 6, units = 'in', res=150)
# here you can pick your TotFish units increase threshold you are interested in (units g per m2) (but results can be intepreted as kg/ha as well)
# where 1 g/m2 = 10 kg/ha
#thres <- c(1.125,1.25,1.5,2,3)
#thres <- c(1.25,1.5,2,3)
#thres <- c(2.5,3,3.5,4,4.5,5)
thres <- c(1,1.5,2.5,3,5,6,8,10)
mainstore <- c()

for(i in 1: length(dats)){
  set <- get(dats[i])
  
  #newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH"=seq(0,30, by=10), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=seq(min(set$SITE_SLOPE_400m_c,na.rm=TRUE),max(set$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25), "ISLAND"=NA, "ECOREGION"=NA, "SITE"=NA,"DIVER"=NA)
  newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH"=seq(0,30, by=10), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=seq(min(set$SITE_SLOPE_400m_c,na.rm=TRUE),max(set$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25), re_formula=NA)
  
  # centre and scale depth
  newdat$DEPTH_c <- (newdat$DEPTH - mean(set$DEPTH))/sd(set$DEPTH)
  
  # order predictions by populations status and depth so that they can be extracted easily below
  newdat <- newdat[order(newdat$POP_STATUS,newdat$DEPTH_c),]
  
  # predict from posterior
  #t3 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE)
  t3 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE,re_formula = NA)
  
  #if(hurd[i]==1){
  #  t33 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE,dpar='hu') 
  #  t3 <- (1-t33)*t3
  #}
  
  # extract samples for unpopulated (works because of ordering step)
  tu <- t3[,1:(dim(t3)[2]/2)]
  
  # extract samples for populated (works because of ordering step)
  tp <- t3[,((dim(t3)[2]/2)+1):dim(t3)[2]]
  
  # calculate differences in samples between depths 
  # i.e. difference in TotFish between 10 and 0m, 20 and 10m, 30 and 20m etc
  
  # create storage 
  tdif <- matrix(NA, nrow=(dim(tu)[2]/length(unique(newdat$DEPTH)))*dim(tu)[1], ncol=(length(unique(newdat$DEPTH))-1))
  tdifp <- matrix(NA, nrow=(dim(tp)[2]/length(unique(newdat$DEPTH)))*dim(tp)[1], ncol=(length(unique(newdat$DEPTH))-1))
  depnms<- c()
  
  # use loop to split samples matrices by depth (works because of order step earlier)
  # each iteration will calculate the difference between two depths 
  for(j in 1:(length(unique(newdat$DEPTH))-1)){
    # unpopulated first
    lower <- as.vector(tu[,(((j-1)*(dim(tu)[2]/length(unique(newdat$DEPTH))))+1):(j*(dim(tu)[2]/length(unique(newdat$DEPTH))))]) # extract all samples at the lower depth
    upper <- as.vector(tu[,((j*(dim(tu)[2]/length(unique(newdat$DEPTH))))+1):((1+j)*(dim(tu)[2]/length(unique(newdat$DEPTH))))]) # extract all samples at the upper depth
    
    #tdif[,j] <- upper/lower   # take difference between samples (samples are still paired, and therefore obtained under the same conditions)
    tdif[,j] <- upper-lower
    
    # and for populated 
    lowerp <- as.vector(tp[,(((j-1)*(dim(tp)[2]/length(unique(newdat$DEPTH))))+1):(j*(dim(tp)[2]/length(unique(newdat$DEPTH))))])
    upperp <- as.vector(tp[,((j*(dim(tp)[2]/length(unique(newdat$DEPTH))))+1):((j+1)*(dim(tp)[2]/length(unique(newdat$DEPTH))))])
    
    #tdifp[,j] <- upperp/lowerp 
    tdifp[,j] <- upperp-lowerp 
    
    # get names of depths being compared for ease of reading later
    depnms <- c(depnms, paste0(unique(newdat$DEPTH)[j], " to ", unique(newdat$DEPTH)[j+1]))
    
  }
  
  tdif <- apply(tdif,2, FUN=sort) # this isn't really necessary but puts the calculated distances in numerical order
  tdifp <- apply(tdifp,2, FUN=sort)
  
  
  # loop through the thresholds to get proportion of differences 
  
  thresstore <- matrix(NA, nrow=length(thres)*2, ncol=dim(tdif)[2])
  
  for(k in 1: length(thres)){
    probz <- apply(tdif,2,FUN=function(x5) {length(x5[which(x5>=thres[k])])/length(x5)}) 
    probzp <- apply(tdifp,2,FUN=function(x5) {length(x5[which(x5>=thres[k])])/length(x5)})
    
    thresstore[((k*2)-1),] <- probz
    thresstore[k*2,] <- probzp
    
  }
  
  thresstore <- as.data.frame(thresstore)
  colnames(thresstore) <- depnms
  
  thresstore$POP_STATUS <- rep(levels(set$POP_STATUS),times=length(thres))
  thresstore$thres <- rep(thres, each=length(levels(set$POP_STATUS)))
  thresstore$model <- rep(dats[i], dim(thresstore)[1])
  
  mainstore <- rbind(mainstore, thresstore)
  
}

#### now plot (Figure 3)

mainstore$POP_STATUS <- factor(mainstore$POP_STATUS, levels=c("U","P"))
mainstore$thres <- factor(mainstore$thres, levels=rev(thres))
ylims <- range(mainstore[,1:(dim(mainstore)[2]-3)])

par(mfrow=c(length(unique(mainstore$POP_STATUS)),length(unique(mainstore$thres))), oma=c(1,1,3,1), mar=c(4,4,2,1))

for(i in 1:length(levels(mainstore$POP_STATUS))){
  set <- mainstore[which(mainstore$POP_STATUS==levels(mainstore$POP_STATUS)[i]),]
  set <- droplevels(set)
  
  for(j in 1:length(levels(mainstore$thres))){
    set2 <- set[which(set$thres==levels(mainstore$thres)[j]),]
    set2 <- droplevels(set2)
    
    #plot(NA, xlim=c(1,dim(mainstore)[2]-3), ylim=c(0,1),las=1, main=paste((as.numeric(as.character(set2$thres[1]))-1)*100,"% increase"), xlab="", ylab="", xaxt='n')
    plot(NA, xlim=c(1,dim(mainstore)[2]-3), ylim=c(0,1),las=1, main=paste(as.numeric(as.character(set2$thres[1]))*10,"kg/ha increase"), xlab="", ylab="", xaxt='n')
    xs <- (1:(dim(set2)[2]-3))
    for(k in 1:dim(set2)[1]){
      ys <- unlist(set2[k,1:(dim(set2)[2]-3)])
      lines(ys~xs, lwd=0.2, col=mycols[k])
      abline(h=0.5, col="gray55",lty=3)
      abline(h=0.75, col="gray55",lty=2)
      abline(h=0.95, col="gray55",lty=1)
      points(ys~xs, pch=16, cex=2, col=mycols[k])
      axis(side=1, at=xs, labels=paste0(names(set2[,1:(dim(mainstore)[2]-3)])))
    }
    
  }
}

mtext("Depth bin (m)", side=1, line=-0.5, outer=TRUE)
mtext("Probability of increase", side=2, line=-0.5, outer=TRUE)
mtext("Unpopulated", side=2, at=0.75, line=-0.5, outer=TRUE)
mtext("Populated", side=2, at=0.25, line=-0.5, outer=TRUE)
#legend.fun("top", legend=nms, pch=16, col=mycols, bty='n', ncol=length(nms), cex=1.5)
dev.off()
