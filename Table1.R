# R FIGURE code for Richardson et al 2022
# Finalized on 26th July 2022

# By LE Richardson, AJ Delargy and P Neubauer 
library(brms)
# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

###################################################################
###################################################################

# Table 1

# this is the same as above but now the threshold is 0, i.e. we want the probability of any increase 

# storage

unpopstore <- matrix(NA, nrow=length(seq(0,30, by=10))-1, ncol=length(dats)) # this one will contain the prob of any increase for unpopulated
popstore <- unpopstore   # # this one will contain the prob of any increase for populated
unpopstore_dec <- unpopstore    # # this one will contain the prob of any decrease for unpopulated
popstore_dec <- unpopstore      # # this one will contain the prob of any decrease for populated

colnames(unpopstore) <- nms
colnames(popstore)<- nms

colnames(unpopstore_dec) <- nms
colnames(popstore_dec)<- nms

# fill em in

for(i in 1: length(dats)){
  set <- get(dats[i])
  
  #newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH"=seq(0,30, by=10), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=seq(min(set$SITE_SLOPE_400m_c,na.rm=TRUE),max(set$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25), "ISLAND"=NA, "ECOREGION"=NA, "SITE"=NA,"DIVER"=NA)
  newdat <- expand.grid("POP_STATUS"=levels(set$POP_STATUS), "DEPTH"=seq(0,30, by=10), "OBS_YEAR"=levels(set$OBS_YEAR), "SITE_SLOPE_400m_c"=seq(min(set$SITE_SLOPE_400m_c,na.rm=TRUE),max(set$SITE_SLOPE_400m_c,na.rm=TRUE),length.out=25), re_formula=NA)
  
  # centre and scale depth
  newdat$DEPTH_c <- (newdat$DEPTH - mean(set$DEPTH))/sd(set$DEPTH)
  
  # order predictions by populations status and depth so that they can be extracted easily below
  newdat <- newdat[order(newdat$POP_STATUS,newdat$DEPTH_c),]
  
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
  
  # calculate differences in samples between depths 
  # i.e. difference in TotFish between 5 and 0m, 10 and 5m, 15 and 10m etc
  
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
    
    tdif[,j] <- upper - lower   # take difference between samples (samples are still paired, and therefore obtained under the same conditions)
    
    # and for populated 
    lowerp <- as.vector(tp[,(((j-1)*(dim(tp)[2]/length(unique(newdat$DEPTH))))+1):(j*(dim(tp)[2]/length(unique(newdat$DEPTH))))])
    upperp <- as.vector(tp[,((j*(dim(tp)[2]/length(unique(newdat$DEPTH))))+1):((j+1)*(dim(tp)[2]/length(unique(newdat$DEPTH))))])
    
    tdifp[,j] <- upperp - lowerp  
    
    # get names of depths being compared for ease of reading later
    depnms <- c(depnms, paste0(unique(newdat$DEPTH)[j], " to ", unique(newdat$DEPTH)[j+1]))
    
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

# could just have calculated the decrease as 1 - the probability of increase: 

unpopstore + unpopstore_dec

popstore + popstore_dec
