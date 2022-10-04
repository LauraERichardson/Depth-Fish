library(tidyverse)
library(brms)
library(tidybayes)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

########## Figure S7
#slope and population status: ggplot version

pphus <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  newdat <- data.frame(expand.grid("POP_STATUS"=levels(set$POP_STATUS), 
                                   "DEPTH_c"=mean(set$DEPTH_c), 
                                   "SITE_SLOPE_400m"=seq(min(set$SITE_SLOPE_400m), max(set$SITE_SLOPE_400m), l=100)))
  newdat$SITE_SLOPE_400m_c <- (newdat$SITE_SLOPE_400m - mean(fish$SITE_SLOPE_400m))/sd(fish$SITE_SLOPE_400m)
  newdat$trophic_group <- factor(nms[i], levels = nms)
  
  newdat %>% 
    add_epred_draws(get(models[i]), 
                    re_formula = NA, seed = 123, dpar=TRUE)
  
}, mc.cores = 5) %>% bind_rows()

pp_bds <- pphus %>% 
  group_by(trophic_group,  POP_STATUS, SITE_SLOPE_400m) %>%
  median_qi(.epred, .width = 0.75)

# data std to 2010
phus <- parallel::mclapply(1:length(models), function(i) {
  set <- get(dats[i])
  pp <- posterior_epred(get(models[i]), 
                        newdata = get(models[i])$data %>% 
                          mutate(DEPTH_c=mean(set$DEPTH_c)))
  
  data.frame(dens = posterior_summary(pp)[,1],
             SITE_SLOPE_400m = get(models[i])$data$SITE_SLOPE_400m*sd(fish$SITE_SLOPE_400m) + mean(fish$SITE_SLOPE_400m),
             POP_STATUS = get(models[i])$data$POP_STATUS,
             trophic_group = factor(nms[i], levels = nms))
  
}, mc.cores = 5) %>% bind_rows()

# restrict axis to plotted intervals
phus <- phus %>% inner_join(pp_bds %>% 
                              group_by(trophic_group) %>% 
                              summarise(md = max(.upper))) %>%
  group_by(trophic_group) %>%
  filter(dens<=md*1.05)

ggplot() + 
  geom_point(aes(x=SITE_SLOPE_400m, y=dens*10,col=trophic_group), alpha=0.2, data=phus%>% filter(POP_STATUS=='U')) +
  geom_point(aes(x=SITE_SLOPE_400m, y=dens*10), col='grey50', alpha=0.05, data=phus%>% filter(POP_STATUS=='P')) +
  geom_ribbon(aes(x=SITE_SLOPE_400m, fill=trophic_group, ymin=.lower*10, ymax=.upper*10),  alpha=0.5, data=pp_bds %>% filter(POP_STATUS=='U')) +
  geom_ribbon(aes(x=SITE_SLOPE_400m, ymin=.lower*10, ymax=.upper*10),  alpha=0.5, fill='grey50', data=pp_bds %>% filter(POP_STATUS=='P')) +
  geom_line(aes(x=SITE_SLOPE_400m, col=trophic_group, y=.epred*10, group=POP_STATUS), data=pp_bds %>% filter(POP_STATUS=='U')) +
  geom_line(aes(x=SITE_SLOPE_400m, y=.epred*10, group=POP_STATUS), col='grey50', data=pp_bds %>% filter(POP_STATUS=='P')) +
  scale_fill_manual('',values = mycols, guide='none') +
  scale_colour_manual('',values = mycols, guide='none') +
  facet_wrap(~trophic_group, scales='free', nrow = 1) +
  xlab("Bathymetric steepness (°)") + 
  ylab("Fish biomass (kg/ha)") +
  cowplot::theme_cowplot() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(face = 'bold')
  )


ggsave('FigureS7_gg.png',width = 10, height = 5, units = 'in',dpi = 150)


########## Figure S7

png('FigureS7.png',width = 7, height = 6, units = 'in', res=150)
#slope and population status:

# name of response variable used in model
responses <- c("TotFish","PRIMARY","PLANKTIVORE","SECONDARY","PISCIVORE")

#split the screen

par(mfrow=c(1,length(models)), oma=c(1,1,3,1), mar=c(4,4,2,1))

for(i in 1:length(dats)){
  set <- get(dats[i])
  
  # remove rows with NA in site-slope 
  set <- set[which(!(is.na(set$SITE_SLOPE_400m_c))),]
  set <- droplevels(set)
  
  #newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH"=seq(min(set$DEPTH),max(set$DEPTH), length.out=25), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=mean(set$SITE_SLOPE_400m_c), "ISLAND"=NA, "ECOREGION"=NA, "SITE"=NA, "DIVER"=NA)
  newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH_c"=mean(set$DEPTH_c), "OBS_YEAR"=levels(set$OBS_YEAR)[1], "SITE_SLOPE_400m"=seq(min(set$SITE_SLOPE_400m),max(set$SITE_SLOPE_400m), length.out=25))
  
  # prepare data for predictions with slope fixed at mean. First two years are used because of interactions involving year 
  
  newdat$SITE_SLOPE_400m_c <- (newdat$SITE_SLOPE_400m - mean(fish$SITE_SLOPE_400m))/sd(fish$SITE_SLOPE_400m)
  # order predictions by populations status and depth so that they can be extracted easily below
  newdat <- newdat[order(newdat$POP_STATUS,newdat$SITE_SLOPE_400m_c),]
  
  # get model info  
  fun_link <- get(models[i])$family$linkfun
  inv_fun <- get(models[i])$family$linkinv
  
  # get residuals 
  #set$res <- residuals(get(models[i]))[,1]
  
  # fitted values 
  t8 <- posterior_epred(get(models[i]), 
                        newdata = set , 
                        re_formula = NULL,
                        allow_new_levels=F)
  
  #set$fitted <- posterior_summary(t8)[,1]
  
  set$partials <- posterior_summary(t8)[,1]
  # copied from ggEffects/residualize_over_grid source code 
  
  # convert from g/m2 to kg/ha for plotting 
  set$partials <- set$partials*10
  
  # subset to unpopulated and populated 
  
  setu <- set[which(set$POP_STATUS=="U"),]
  setp <- set[which(set$POP_STATUS=="P"),]
  setu <- droplevels(setu)
  setp <- droplevels(setp)
  
  # to make a nice plot, need to go back to response predictions and pool across other variables (i.e. year and slope)
  
  ##here we need to add the posterior epred of the probability of observing a zero
  ##The expected value for total model is (1-E[p])*E[d]
  
  # predict from posterior
  #t3  <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE) 
  t3  <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE, re_formula = NA)
  
  #if(hurd[i]==1){
  #  t33 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE,dpar='hu') 
  #  t3 <- (1-t33)*t3
  #}
  
  # convert to kg/ha
  t3 <- t3*10 
  
  # extract samples for unpopulated (works because of ordering step)
  tu <- t3[,1:(dim(t3)[2]/2)]
  
  # extract samples for populated (works because of ordering step)
  tp <- t3[,((dim(t3)[2]/2)+1):dim(t3)[2]]
  
  depthsu <- matrix(NA, nrow=length(unique(newdat$SITE_SLOPE_400m)), ncol=8)
  depthsp <- matrix(NA, nrow=length(unique(newdat$SITE_SLOPE_400m)), ncol=8)
  
  
  for(j in 1:length(unique(newdat$SITE_SLOPE_400m))){ # loop through each depth and calculate median and CIs
    loweru <- as.vector(tu[,(((j-1)*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):(j*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))]) # extract all samples at specified depth at unpopulated sites
    lowerp <- as.vector(tp[,(((j-1)*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):(j*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))]) # extract all samples at specified depth at populated sites
    
    qu <- quantile(loweru, probs=c(0.025, 0.125,0.25, 0.75, 0.875,0.975))
    qp <- quantile(lowerp, probs=c(0.025, 0.125,0.25, 0.75,0.875,0.975))
    
    depthsu[j,1] <- unique(newdat$SITE_SLOPE_400m)[j]
    depthsu[j,2:7] <- qu
    depthsu[j,8] <- median(loweru)
    
    depthsp[j,1] <- unique(newdat$SITE_SLOPE_400m)[j]
    depthsp[j,2:7] <- qp
    depthsp[j,8] <- median(lowerp)
    
  }
  
  # plot
  depthst <- rbind(depthsu,depthsp)# for plot size
  
  #yvals <- depthst[,2:dim(depthst)[2]]  # 95% CIs 
  yvals <- depthst[,3:(dim(depthst)[2]-2)]  # 75% CIs
  #yavls <- depthst[,4:(dim(depthst)[2]-3)] # 50% CIs
  
  #yvals <- c(stack(as.data.frame(depthst[,2:dim(depthst)[2]]))[,1],setu$partials,setp$partials)
  
  # THE ABOVE LINES CONTROL THE Y LIMITS ON EACH PANEL (yvals)
  # the first three lines will curtail the upper y limit at the confidence intervals, and therefore some partial residuals will be excluded 
  # the second line will include a range large enough to get all partial residuals 
  # to use one, hash out the others 
  
  plot(depthsp[,dim(depthsp)[2]]~depthsp[,1], type='l', lwd=0, las=1, xlab="", ylab="", ylim=c(0,max(yvals)), main=nms[i])
  points(setu$partials~setu$SITE_SLOPE_400m, col=adjustcolor(mycols[i], alpha.f=0.5), pch=16)
  points(setp$partials~setp$SITE_SLOPE_400m, col=adjustcolor("gray55", alpha.f=0.5), pch=16)
  
  # IN THE BELOW, USE THE alpha.f argument to set the transparency (or degree of shading), where 0 is entirely transparent and 1 is not
  
  #polygon(c(depthsu[,1],rev(depthsu[,1])),c(depthsu[,2],rev(depthsu[,7])),border=NA,col=adjustcolor(mycols[i], alpha.f=0.2))# 95% CIs
  polygon(c(depthsu[,1],rev(depthsu[,1])),c(depthsu[,3],rev(depthsu[,6])),border=NA,col=adjustcolor(mycols[i], alpha.f=0.2))# 75% CIs
  #polygon(c(depthsu[,1],rev(depthsu[,1])),c(depthsu[,4],rev(depthsu[,5])),border=NA,col=adjustcolor(mycols[i], alpha.f=0.2))# 50% CIs
  lines(depthsu[,dim(depthsu)[2]]~depthsu[,1], lwd=3, col=mycols[i])
  
  #polygon(c(depthsp[,1],rev(depthsp[,1])),c(depthsp[,2],rev(depthsp[,7])),border=NA,col=adjustcolor("gray55", alpha.f=0.5))# 95% CIs
  polygon(c(depthsp[,1],rev(depthsp[,1])),c(depthsp[,3],rev(depthsp[,6])),border=NA,col=adjustcolor("gray55", alpha.f=0.5))# 75% CIs
  #polygon(c(depthsp[,1],rev(depthsp[,1])),c(depthsp[,4],rev(depthsp[,5])),border=NA,col=adjustcolor("gray55", alpha.f=0.5))# 50% CIs
  lines(depthsp[,dim(depthsp)[2]]~depthsp[,1], lwd=3, col="gray48")
  
}

mtext("Bathymetric steepness (°)", side=1, line=-1, outer=TRUE)
mtext("Fish biomass (kg/ha)", side=2,line=-1, outer=TRUE)
#legend.fun("top", legend=c("Populated", rep("Unpopulated", length(models))), lty=1, lwd=3, col=c("gray55",mycols), bty='n', cex=1.2, ncol=length(models)+1)
dev.off()
####
