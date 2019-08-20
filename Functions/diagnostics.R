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
  
  shp1$predictions  = glmpreds_tmp[match( shp1$GID_1, dat$adm1)]
  
  
  mybreaks = seq(0, 1.0001, length.out=1001)
  mycols =  colours
  vcols = findInterval(shp1$predictions, mybreaks)
  
  
  plot(shp0)
  plot(shp1, col=mycols[vcols] , lty=0, add=TRUE)
  plot(shp0, add=TRUE)
  fields::image.plot(legend.only=TRUE, breaks=mybreaks, col=mycols, zlim=c(0,1), horizontal = TRUE)
  
}

#------------------------------------------------------------------------------------------------------
plot_sero = function(seroout,
                     hpd_sero,
                     pop_agg3d,
                     dim_year,
                     dim_age,
                     inc_v3d,
                     pop_moments_agg,
                     vc_agg3d){
  
  
  
  p = create_pop_at_survey(pop_agg3d,
                           seroout$sero_studies,
                           dim_year)
  
  p_at_survey_3d=p$p_at_survey_3d
  P_tot_survey_2d=p$P_tot_survey_2d
  
  
  for (i in 1:3){
    
    #set parameters
    param =  hpd_sero[,i]
    v_e = hpd_sero[1, i]
    
    seroprev_predict_survey_r = seroprev_predict_survey_f =  NULL
    
    for (index_survey in 1:seroout$no_sero_surveys){
      
      vcfac = ifelse(!is.na(seroout$vc_factor[index_survey]),
                     seroout$vc_factor[index_survey],
                     param[grep("vc_fac", names( param))])
      
      # SEROPREVALENCE FOI #
      Foisurv = param[grep("^Foi", names( param))][index_survey]
      
      vc_agg_ag = vc_agg3d[index_survey, paste0(seroout$study_years[index_survey]),]
      vc=vcfac*v_e*vc_agg_ag
      
      pop = pop_agg3d[index_survey,paste0(seroout$study_years[index_survey]),]
      
      sero1 = c(as.matrix((1-exp(-Foisurv*0:100))*pop))
      sero_ag = aggregate(sero1, by=list(age_group=0:100), sum)
      sero_ag$x = sero_ag$x/aggregate(pop, by=list(age_group=0:100), sum)$x
      sero_ag$x[is.na(sero_ag$x)]=0
      
      res_out = 1 - (1-sero_ag$x[sero_ag$age_group>0])*(1-vc)
      
      seroprev_predict_survey_f = rbind(seroprev_predict_survey_f,  res_out)
      
    } # end of for(index_survey)
    seroprev_predict_survey_print_f = data.frame(sero_studies=rep(seroout$sero_studies, times = 1),
                                                 seroprev_predict_survey_f)
    
    
    if (i==1) {
      f_sero_predictions_lo=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }
    if (i==2) {
      f_sero_predictions=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }
    if (i==3) {
      f_sero_predictions_hi=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }
    
    
  } #end of for(i)
  
  
  ### PLOT ###
  par(mfrow=c(6,7),  oma = c(4, 4, 0, 0), mar = c(2, 2, 1, 1))
  for (i in 1:seroout$no_sero_surveys){
    #get upper bound for plot axes
    maxsero = max(as.numeric(f_sero_predictions_hi[i,1:86]),
                  na.rm = TRUE)
    
    #pull out key parts of data for plotting
    age_points = ( seroout$survey_dat[[i]]$age_min+seroout$survey_dat[[i]]$age_max )/2
    seroprevalence = seroout$survey_dat[[i]]$positives / seroout$survey_dat[[i]]$samples
    
    #plot the data
    plot( age_points,
          seroprevalence,
          main = paste0((seroout$sero_studies[i])),
          cex = 1,
          ylab = "Seroprevalence", pch=19, col = alpha("black",0.5),
          xlim=c(0,85), ylim=c(0,maxsero+0.1 ), xlab = "Age")
    
    #add the HPD of predictions
    polygon(c(85:0,0:85), c(rev(as.numeric(f_sero_predictions_hi[i,1:86])),
                            (as.numeric(f_sero_predictions_lo[i,1:86]))) ,
            border=NA, col = rgb(0, 0, 1,0.2) )
    
    #add the median predictions
    lines(0:85,as.numeric(f_sero_predictions[i,1:86]),
          type="l", main = paste0((seroout$sero_studies[i])),
          ylab = "Seroprevalence", col = "blue" )
    
    ## add binomial confidence for the data
    conf_int = Hmisc::binconf(seroout$survey_dat[[i]]$positives,
                              seroout$survey_dat[[i]]$samples)
    
    arrows(age_points,
           conf_int[,2],
           age_points,
           conf_int[,3], length=0.05, angle=90, code=3)
    
    points(age_points, conf_int[,1])
  }
  mtext('Age', side = 1, outer = TRUE, line = 2)
  mtext('Seroprevalence', side = 2, outer = TRUE, line = 2)
}

#------------------------------------------------------------------------------------------------------#
get_hpd = function(mcmc_out,
                   glm_mcmc_out){
  
  if(!anyNA(mcmc_out)){
    ### hpd for serology ###
    #full uncertainty for shared parameters
    hpd_sero = hpd_out(mcmc_out[,grep("^vac", names(mcmc_out))])
    
    
    #hpd by model type
    hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out[,grep("^Foi", names(mcmc_out))]) )
    
    #last shared parameter
    hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out[,grep("^vc", names(mcmc_out))]))
  }
  
  if(!anyNA(glm_mcmc_out)){
    ### hpd for glm ###
    hpd_glm = hpd_out( subset( glm_mcmc_out, select = -c(posteriorProb, acceptRate)),
                       log_scale = FALSE)
  }
  
  if(!anyNA(mcmc_out) & !anyNA(glm_mcmc_out)){
    
    rownames(hpd_sero)[c(1, nrow(hpd_sero))] = c("vac_eff", "vc_factor_CMRs")
    
    ### bind all together for each model ###
    R0_param = rbind(hpd_sero[grep("^vac", rownames(hpd_sero)),],
                     hpd_glm,
                     hpd_sero[ grep("^R0", rownames(hpd_sero) ),],
                     hpd_sero[grep("^vc", rownames(hpd_sero)),])
    
    Foi_param = rbind(hpd_sero[grep("^vac", rownames(hpd_sero)),],
                      hpd_glm,
                      hpd_sero[ grep("^Foi", rownames(hpd_sero) ),],
                      hpd_sero[grep("^vc", rownames(hpd_sero)),])
    
    rownames(R0_param)[c(1,nrow(R0_param))] =
      rownames(Foi_param)[c(1,nrow(Foi_param))] = c("vac_eff", "vc_factor_CMRs")
    
    
    
    out = list(Foi_param = Foi_param, R0_param = R0_param, hpd_sero = hpd_sero, hpd_glm)
  } else if(anyNA(glm_mcmc_out)){
    rownames(hpd_sero)[c(1, nrow(hpd_sero))] = c("vac_eff", "vc_factor_CMRs")
    out = list(hpd_sero = hpd_sero)
  } else {
    out = list(hpd_glm = hpd_glm)
  }
  
  return(out)
}

#------------------------------------------------------------------------------------------------------#
plot_transmission_intensity2 = function(x,
                                        ii,
                                        seroout,
                                        params,
                                        dat,
                                        t0_vac_africa,
                                        dim_year,
                                        dim_age,
                                        p_prop_3d,
                                        P_tot_2d,
                                        inc_v3d,
                                        pop1,
                                        vc2d,
                                        varsin_nc,
                                        polydeg,
                                        R0_lookup,
                                        model_type,
                                        shp1,
                                        shp0,
                                        colours){
  
  
  runs = fun_calc_transmission_Africa(x,
                                      ii ,
                                      seroout ,
                                      params ,
                                      dat ,
                                      t0_vac_africa ,
                                      dim_year ,
                                      dim_age ,
                                      p_prop_3d ,
                                      P_tot_2d ,
                                      inc_v3d ,
                                      pop1,
                                      vc2d,
                                      varsin_nc,
                                      polydeg,
                                      R0_lookup,
                                      model_type)
  
  runs = switch(model_type,
                "R0" = runs -1,
                "Foi" = runs)
  
  shp1$runs = NA 
  shp1$runs = runs[match(shp1$GID_1,dat$adm0_adm1)]
  
  mybreaks= seq(min(log10(shp1$runs), na.rm = TRUE), max(log10(shp1$runs), na.rm = TRUE)+0.01, length.out=101)

  vcols = findInterval(log10(shp1$runs),mybreaks)
  
  plot(shp0)
  plot(shp1,col=colours[vcols], lty=0, add=T)
  
  plot(shp0, lwd=2, add=T)
  
  
  
  image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
             axis.args = list(at = c(-5:1), labels =c("1e-05","1e-04", "0.001", "0.01", "0.1","1","10"), las =2),
             legend.mar = 3.5)
  
}
