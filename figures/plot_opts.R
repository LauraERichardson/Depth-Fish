# global legend function

legend.fun<-function(...){
  opar<-par(fig=c(0,1,0,1),oma=rep(0,4),
            mar=rep(0,4), new=TRUE)
  on.exit(par(opar))
  plot(0,0,type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

# plot display names 
nms <- c("Total biomass", "Primary consumer","Planktivore","Secondary consumer","Piscivore")
# colour for each model
mycols <- c("#0072B2","#009E73","#56B4E9","#CC79A7","#E69F00")
hurd  <- c(F,F,T,F,T)
#desired_order <- c("Depth","Slope","Pop_Status","Depth:Pop_Status","Depth Slope smooth 1","Depth Slope smooth 2") 
# in the above line you can specify the order of fixed effects on the y-axis of the below plot 
# names of of each variable must exactly match those used to label the row names of q3 below 
