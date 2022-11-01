
# By LE Richardson, AJ Delargy and P Neubauer 

require(tidyverse)
library(brms)
if(!require(Ternary)) install.packages('Ternary')
require(Ternary)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############## display percentage of random effect variance from each other as a ternary plot 

pdf('Figure5.pdf',width = 7, height = 6)

# Figure 5
gmean <- function(x) exp(mean(log(x)))

c.pan <- floor(sqrt(length(models)))	# automatically determine number of cols to split screen
r.pan <- ceiling(length(models)/c.pan)	# automatically determine rows
par(mfrow=c(r.pan,c.pan), mar=c(1,1,1,1))

for(j in 1:length(models)){
  v1 <- posterior_samples(get(paste(models[j])), pars='sd')
  
  per5 <- data.frame(Ecoregion = v1$sd_ECOREGION__Intercept, Island = v1$sd_ISLAND__Intercept, Site = v1$sd_SITE__Intercept)
  
  if(hurd[j]==1){
    per5 <- data.frame(Ecoregion = v1$sd_ECOREGION__Intercept+v1$sd_ECOREGION__hu_Intercept, Island = v1$sd_ISLAND__Intercept + v1$sd_ISLAND__hu_Intercept, Site = v1$sd_SITE__Intercept + v1$sd_SITE__hu_Intercept)
  }
  
  # plot 
  TernaryPlot(alab = "ECOREGION", blab = "ISLAND", clab = "SITE",
              point = 'up', 
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate =FALSE,
              padding = 0.08, clockwise=TRUE, lab.cex=1.25, axis.cex=1.25, main=nms[j])
  
  colfunc <- colorRampPalette(c("white",paste(mycols)[j]))
  yes <- colfunc(40)
  ColourTernary(TernaryDensity(per5, resolution = 40L),spectrum = adjustcolor(yes, alpha.f=0.7))
  TernaryPoints(apply(per5,2,FUN=gmean), col='black', pch=16)
  dum1 <- apply(per5,2,FUN=gmean)
  dum2 <- dum1/sum(dum1)
  TernaryArrows(dum1, c(1-dum2[2],dum2[2],0), length = 0.02, col ='black')
  TernaryArrows(dum1, c(dum2[1],0,1-dum2[1]), length = 0.02, col ='black')
  TernaryArrows(dum1, c(0,1-dum2[3],dum2[3]), length = 0.02, col ='black')
  #strip.background = element_rect(colour = "black", fill = "white"),
  #strip.text.x = element_text(face = 'bold')
}

dev.off()


# Table S10

q5 <- c()

for(j in 1:length(models)){
  v1 <- posterior_samples(get(paste(models[j])), pars='sd')
  
  per5 <- data.frame(Ecoregion = v1$sd_ECOREGION__Intercept, Island = v1$sd_ISLAND__Intercept, Site = v1$sd_SITE__Intercept)
  
  if(hurd[j]==1){
    per5 <- data.frame(Ecoregion = v1$sd_ECOREGION__Intercept+v1$sd_ECOREGION__hu_Intercept, Island = v1$sd_ISLAND__Intercept + v1$sd_ISLAND__hu_Intercept, Site = v1$sd_SITE__Intercept + v1$sd_SITE__hu_Intercept)
  }
  
  # get medians and CIs 
  
  q4 <- t(apply(per5, 2, FUN=function(x1) {quantile(x1, probs=c(0.5, 0.025, 0.125, 0.875, 0.975))}))
  
    # store 
  q4 <- as.data.frame(q4)
  q4$Model <- rep(nms[j],dim(q4)[1])
  q4$Spatial <- rownames(q4)
  q4 <- q4[order(rev(q4$Spatial)),]
  q5 <- rbind(q5,q4)
  
}

write_csv(as.data.frame(q5), 'tableS10.csv')
