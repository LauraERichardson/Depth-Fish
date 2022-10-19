# R FIGURE code for Richardson et al 2022
# Finalized on 26th July 2022

# By LE Richardson, AJ Delargy and P Neubauer 

library(tidyverse)
library(brms)

# load models and data from common script
source('load_data_and_models_for_figs.R')
source('plot_opts.R')

# get fixed effects samples from posterior and calculate CIs (as quantiles)
png('Figure1.png',width = 7, height = 6, units = 'in', res=150)
allcis <- c()

for(i in 1:length(models)){
  v1 <- posterior_samples(get(paste(models[i])))
  
  # these are scaled because they were fit as scaled in the model 
  q3 <- t(apply(v1,2,FUN=function(x1) {quantile(x1, probs=c(0.025, 0.125, 0.875, 0.975, 0.5))}))  # CIs as quantiles
  
  # subset to fixed effects
  q3 <- q3[which(substr(rownames(q3),1,1)=="b"),]
  
  # cut out year 
  q3 <- q3[which(!(substr(rownames(q3),3,5)=="OBS")),]
  q3 <- q3[which(!(substr(rownames(q3),6,8)=="OBS")),]
  
  # cut out intercept term(s) 
  q3 <- q3[which(!(substr(rownames(q3),3,5)=="Int")),]
  q3 <- q3[which(!(substr(rownames(q3),6,8)=="Int")),]
  
  
  q3 <- as.data.frame(q3)
  # save variable names 
  q3$var <- rownames(q3)
  q3$troph <- rep(substr(models[i],1,3),dim(q3)[1])
  allcis <- rbind(allcis, q3)
}

# reorder for plotting 
allcis$troph <- factor(allcis$troph, levels=unique(allcis$troph))

# remove smooths
allcis <- allcis[-grep('bs',allcis$var),]

# rename the variables 
allcis$var <- as.factor(allcis$var)
levels(allcis$var) <- c("Depth",
                        "Depth:Pop Status", 
                        "Depth hurdle",
                        "Depth:Pop Status hurdle",
                        "Pop Status hurdle",
                        "Steepness hurdle",
                        "Pop Status",
                        "Steepness")


desired_order <- c("Depth:Pop Status hurdle","Pop Status hurdle","Steepness hurdle","Depth hurdle","Depth:Pop Status","Pop Status","Steepness","Depth") 
allcis$var <- factor(allcis$var, levels=desired_order)
allcis$var_num <- as.numeric(allcis$var)
allcis <- allcis[order(allcis$troph,allcis$var),]

### PLOT

lineadj <- rev(seq(-0.35,0.35, length.out=length(models)))  # here you can control the gap between points
# keep the two numbers identical (other than sign) to create evenly spaced groups
# smaller numbers which create closer points and large will create more spread out points 

par(mar=c(5,15,1,1))# this controls the size of margins required to show the labels
par(mfrow=c(1,1))

plot(NA, pch=16, xlab="Scaled fixed effects estimates",  yaxt='n', ylab='', 
     xlim=range(allcis[,1:5]), 
     ylim=range(1:length(levels(allcis$var))) + range(lineadj))
abline(v=0, lty=2, lwd=2, col="grey")       # add line to highlight 0 

for(i in 1:length(unique(allcis$troph))){
  set <-   allcis[which(allcis$troph==levels(allcis$troph)[i]),]
  set <- droplevels(set)
  points((set$var_num + lineadj[i]) ~ set[,5], pch=16, col=mycols[i], cex=2)
  arrows(x0=set[,2], y0=set$var_num+lineadj[i], x1= set[,3], y1=set$var_num+lineadj[i], length=0, lwd=5, col=mycols[i])  # 75% credible intervals    
  arrows(x0=set[,1], y0=set$var_num+lineadj[i], x1= set[,4], y1=set$var_num+lineadj[i], length=0, lwd=2, col=mycols[i])  # 95% credible intervals
}
axis(side=2,at=1:length(desired_order), labels=desired_order, las=1)  # add labels to y-axis
#axis(side=2,at=1:length(desired_order), labels=desired_order, las=1)  # add labels to y-axis
#legend("topleft", legend=nms, pch=16, col=mycols, pt.cex=2, cex=1.5)

dev.off()
