#diagnostics functions

plot_glm_map2 = function(shp0,
                         shp1,
                         dat,
                         Est_beta,
                         x,
                         colours,
                         plot_data = TRUE){
  if(plot_data){
    par(mfrow=c(2,1))
    
    ### data ###
    plot(shp0)
    
    pres= dat$adm1[dat$cases_or_outbreaks>0]
    
    mm1<-match(pres, shp1$GID_1)
    
    
    plot(shp1[mm1,], col=colours[750], add=TRUE)
    plot(shp0,lwd=2, add=TRUE)
    plot(shp1,lwd=1, add=TRUE)
  }
  
  ### model ###
  glmpreds_tmp = fun_calcPred( Est_beta,x,type="response")
  
  mybreaks = seq(0, 1.0001, length.out=1001)
  mycols =  colours
  mm = match(shp1$GID_1, dat$adm1)
  vcols = findInterval(glmpreds_tmp, mybreaks)
  
  
  plot(shp0)
  plot(shp1[!is.na(mm),], col=mycols[vcols] , lty=0, add=TRUE)
  plot(shp0, add=TRUE)
  fields::image.plot(legend.only=TRUE, breaks=mybreaks, col=mycols, zlim=c(0,1), horizontal = TRUE)
  
}
