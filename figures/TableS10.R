# R code for Richardson et al 2022
# Finalized on 26th July 2022

# By LE Richardson, AJ Delargy and P Neubauer 

require(tidyverse)
library(brms)
# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

###################################################################
###################################################################

# Table S10
# Probability of change in biomass with steepness

# storage

unpopstore <- matrix(NA, nrow=length(seq(0,40, by=10))-1, ncol=length(dats)) # this one will contain the prob of any increase for unpopulated
popstore <- unpopstore   # # this one will contain the prob of any increase for populated
unpopstore_dec <- unpopstore    # # this one will contain the prob of any decrease for unpopulated
popstore_dec <- unpopstore      # # this one will contain the prob of any decrease for populated

colnames(unpopstore) <- nms
colnames(popstore)<- nms

colnames(unpopstore_dec) <- nms
colnames(popstore_dec)<- nms


for(i in 1: length(dats)){
  set <- get(dats[i])
  
  newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "SITE_SLOPE_400m"=seq(0,40, by=10), "DEPTH_c"=seq(min(set$DEPTH_c,na.rm=TRUE),max(set$DEPTH_c,na.rm=TRUE),length.out=25), re_formula=NA)
  
  # centre and scale SITE_SLOPE_400m
  newdat$SITE_SLOPE_400m_c <- (newdat$SITE_SLOPE_400m - mean(set$SITE_SLOPE_400m))/sd(set$SITE_SLOPE_400m)
  
  # order predictions by populations status and SITE_SLOPE_400m so that they can be extracted easily below
  newdat <- newdat[order(newdat$POP_STATUS,newdat$SITE_SLOPE_400m_c),]
  
  # predict from posterior
  t3 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE, re_formula = NA)
  
  #if(hurd[i]==1){
  #  t33 <- posterior_epred(get(models[i]), newdata=newdat, allow_new_levels=TRUE,dpar='hu') 
  #  t3 <- (1-t33)*t3
  #}
  
  # extract samples for unpopulated (works because of ordering step)
  tu <- t3[,1:(dim(t3)[2]/2)]
  
  # extract samples for populated (works because of ordering step)
  tp <- t3[,((dim(t3)[2]/2)+1):dim(t3)[2]]
  
  # calculate differences in samples between SITE_SLOPE_400m 
  # i.e. difference in TotFish between 10 and 0 degrees, 20 and 10 degrees, 30 and 20 degrees etc
  
  # create storage 
  tdif <- matrix(NA, nrow=(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m)))*dim(tu)[1], ncol=(length(unique(newdat$SITE_SLOPE_400m))-1))
  tdifp <- matrix(NA, nrow=(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m)))*dim(tp)[1], ncol=(length(unique(newdat$SITE_SLOPE_400m))-1))
  depnms<- c()
  
  # use loop to split samples matrices by SITE_SLOPE_400m (works because of order step earlier)
  # each iteration will calculate the difference between two SITE_SLOPE_400m levels 
  for(j in 1:(length(unique(newdat$SITE_SLOPE_400m))-1)){
    # unpopulated first
    lower <- as.vector(tu[,(((j-1)*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):(j*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))]) # extract all samples at the lower SITE_SLOPE_400m
    upper <- as.vector(tu[,((j*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):((1+j)*(dim(tu)[2]/length(unique(newdat$SITE_SLOPE_400m))))]) # extract all samples at the upper SITE_SLOPE_400m
    
    tdif[,j] <- upper - lower   # take difference between samples (samples are still paired, and therefore obtained under the same conditions)
    
    # and for populated 
    lowerp <- as.vector(tp[,(((j-1)*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):(j*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))])
    upperp <- as.vector(tp[,((j*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))+1):((j+1)*(dim(tp)[2]/length(unique(newdat$SITE_SLOPE_400m))))])
    
    tdifp[,j] <- upperp - lowerp  
    
    # get names of SITE_SLOPE_400m being compared for ease of reading later
    depnms <- c(depnms, paste0(unique(newdat$SITE_SLOPE_400m)[j], " to ", unique(newdat$SITE_SLOPE_400m)[j+1]))
    
  }
  
  tdif <- apply(tdif,2, FUN=sort) # this isn't really necessary but puts the calculated distances in numerical order
  tdifp <- apply(tdifp,2, FUN=sort)
  
  
  thres <- 0
  
  # increase
  probz <- apply(tdif,2,FUN=function(x5) {length(x5[which(x5>thres)])/length(x5)}) 
  probzp <- apply(tdifp,2,FUN=function(x5) {length(x5[which(x5>thres)])/length(x5)})
  
  # decrease 
  probz_dec <- apply(tdif,2,FUN=function(x5) {length(x5[which(x5<thres)])/length(x5)}) 
  probzp_dec <- apply(tdifp,2,FUN=function(x5) {length(x5[which(x5<thres)])/length(x5)})
  
  # store for each model 
  
  unpopstore[,i] <- probz
  popstore[,i] <- probzp
  
  unpopstore_dec[,i] <- probz_dec
  popstore_dec[,i] <- probzp_dec
  
}




rownames(unpopstore) <- depnms
rownames(popstore) <- depnms
rownames(unpopstore_dec) <- depnms
rownames(popstore_dec) <- depnms

# results are in those matrices 
unpopstore
popstore

write_csv(as.data.frame(popstore), 'tableS10.csv')
write_csv(as.data.frame(unpopstore), 'tableS10b.csv')


# could just have calculated the decrease as 1 - the probability of increase: 

unpopstore + unpopstore_dec

popstore + popstore_dec