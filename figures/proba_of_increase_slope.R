library(tidyverse)
library(ggridges)
library(brms)
library(tidybayes)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

############ # Figure 4B Probability of proportional increase  

png('Figure4B.png',width = 3, height = 7, units = 'in', res=150)
# calculate the changes 

thres <- c(1.25,1.5,2)   # 25%, 50%, 100% increase
mainstore <- c()

for(i in 1: length(dats)){
  set <- get(dats[i])
  
  #newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH_c"=seq(min(set$DEPTH_c),max(set$DEPTH_c), length.out=25), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=seq(min(fish$SITE_SLOPE_400m_c,na.rm=TRUE),max(fish$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25), "ISLAND"=NA, "ECOREGION"=NA, "SITE"=NA,"DIVER"=NA)
  newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH_c"=seq(min(set$DEPTH_c),max(set$DEPTH_c), length.out=25),  "SITE_SLOPE_400m_c"=seq(min(fish$SITE_SLOPE_400m_c,na.rm=TRUE),max(fish$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25))
  
  # get slope on natural scale for later 
  newdat$slope <- (newdat$SITE_SLOPE_400m_c*sd(fish$SITE_SLOPE_400m,na.rm=TRUE)) + mean(fish$SITE_SLOPE_400m,na.rm=TRUE)
  
  # order predictions by populations status and depth so that they can be extracted easily below
  newdat <- newdat[order(newdat$SITE_SLOPE_400m_c),]
  
  # predict from posterior
  #t3 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE)
  t3 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE, re_formula=NA)
  
  #if(hurd[i]==1){
  #  t33 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE,dpar='hu') 
  #  t3 <- (1-t33)*t3
  #}
  
  # calculate differences in samples between depths 
  # i.e. difference in TotFish between 5 and 0m, 10 and 5m, 15 and 10m etc
  
  # create storage 
  tdif <- matrix(NA, nrow=(dim(t3)[2]/length(unique(newdat$SITE_SLOPE_400m_c)))*dim(t3)[1], ncol=(length(unique(newdat$SITE_SLOPE_400m_c))-1))
  depnms<- c()
  
  # use loop to split samples matrices by depth (works because of order step earlier)
  # each iteration will calculate the difference between two depths
  
  # extract estimates at initial slopes (lowest values)
  initials <- as.vector(t3[,1:((dim(t3)[2]/length(unique(newdat$SITE_SLOPE_400m_c))))]) # extract all samples at the lower depth
  
  # loop through remaining slope values and calculate ratio between them and initials
  for(j in 1:(length(unique(newdat$SITE_SLOPE_400m_c))-1)){
    
    upper <- as.vector(t3[,((j*(dim(t3)[2]/length(unique(newdat$SITE_SLOPE_400m_c))))+1):((1+j)*(dim(t3)[2]/length(unique(newdat$SITE_SLOPE_400m_c))))]) # extract all samples at the upper depth
    
    tdif[,j] <- upper/initials   # take ratio of samples (samples are still paired, and therefore obtained under the same conditions)
    
  }
  
  # loop through the thresholds to get proportion of differences 
  
  thresstore <- matrix(NA, nrow=length(thres), ncol=dim(tdif)[2])
  
  for(k in 1: length(thres)){
    probz <- apply(tdif,2,FUN=function(x5) {length(x5[which(x5>=thres[k])])/length(x5)}) 
    
    thresstore[k,] <- probz
    
  }
  
  thresstore <- as.data.frame(thresstore)
  
  thresstore$thres <- thres
  thresstore$model <- rep(dats[i], dim(thresstore)[1])
  
  mainstore <- rbind(mainstore, thresstore)
  
}

#### now plot (Figure 4B)

# split screen 

par(mfrow=c(length(nms),1), oma=c(3,3,3,1), mar=c(2,4,2,1))

legendpos <- c("topright","topright","bottom","topright","topright")

shadecols <- matrix(NA, nrow=length(models), ncol=length(thres))

colfunc <- colorRampPalette(c("lightblue1","darkblue"))
shadecols[1,] <-  colfunc(3)
colfunc <- colorRampPalette(c("lightgreen","darkgreen"))
shadecols[2,] <- colfunc(3)
colfunc <- colorRampPalette(c("lightblue1","deepskyblue4"))
shadecols[3,] <- colfunc(3)
colfunc <- colorRampPalette(c("lightpink","magenta4"))
shadecols[4,] <- colfunc(3)
colfunc <- colorRampPalette(c("lightsalmon","darkorange4"))
shadecols[5,] <- colfunc(3)

loopcount <- 1

slopvals <- seq(min(fish$SITE_SLOPE_400m, na.rm=TRUE),max(fish$SITE_SLOPE_400m,na.rm=TRUE),length.out=25)[-1]

for(i in 1:length(unique(mainstore$model))){
  set <- mainstore[which(mainstore$model==unique(mainstore$model)[i]),]
  set <- droplevels(set)
  
  plot(NA, xlim=range(slopvals), ylim=c(0,1),las=1, xlab="", ylab="", xaxt='n')
  xs <- slopvals
  for(k in 1:dim(set)[1]){
    ys <- unlist(set[k,1:(dim(set)[2]-2)])
    abline(h=0.75, col="gray55",lty=2)
    lines(ys~xs, lwd=2, col=shadecols[i,k])
  }  
  axis(side=1, at=seq(floor(slopvals[1]/5)*5,ceiling(slopvals[length(slopvals)]/5)*5,by=5), labels=seq(floor(slopvals[1]/5)*5,ceiling(slopvals[length(slopvals)]/5)*5,by=5))
  #legend(legendpos[loopcount], legend=paste((thres-1)*100,"%"),lty=1, lwd=2, col=shadecols[i,], bty='n')
  loopcount <- loopcount + 1
  
}

#mtext("Bathymetric steepness (degrees)", side=1, line=1, outer=TRUE)
mtext("Probability of increased biomass", side=2, line=1, outer=TRUE)
#mtext(rev(nms), side=2, line=-0.5, at=((1:5)/5)-0.1, outer=TRUE)
dev.off()
###################################################################
###################################################################
