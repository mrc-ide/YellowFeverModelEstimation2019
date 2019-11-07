#diagnostics functions

plot_glm_map2 = function(shp0,
                         shp1,
                         dat,
                         Est_beta,
                         x,
                         colours,
                         plot_data = TRUE){
  if(plot_data){
    par(mfrow=c(2,1), omi=c(0.5,0.3,0,0), plt=c(0.1,0.9,0,0.7))
    
    ### data ###
    plot(shp0)
    
    pres= dat$adm1[dat$cases_or_outbreaks>0]
    
    mm1<-match(pres, shp1$GID_1)
    
    
    plot(shp1[mm1,], col=colours[length(colours)], add=TRUE)
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
      
      res_out = 1 - (1-sero_ag$x[sero_ag$age_group>=0])*(1-vc)
      
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
plot_sero_R0 = function(seroout,
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
    
    seroprev_predict_survey_r = NULL
    
    for (index_survey in 1:seroout$no_sero_surveys){
      
      vcfac = ifelse(!is.na(seroout$vc_factor[index_survey]),
                     seroout$vc_factor[index_survey],
                     param[grep("vc_fac", names( param))])
      
      # SEROPREVALENCE R0 #
      R0surv = param[grep("R0", names( param))][index_survey]
      res = R0_recurrence_seroprev_survey(R0=R0surv,
                                          foi_const = 0,
                                          t0_vac=seroout$t0_vac[index_survey],
                                          dim_year=dim_year,
                                          dim_age=dim_age,
                                          p_at_survey= p_at_survey_3d[index_survey,,],
                                          P_tot_survey= P_tot_survey_2d[index_survey,],
                                          inc_v3d=inc_v3d_agg[index_survey,,],
                                          pop_moments=pop_moments_agg[index_survey,],
                                          vac_eff_arg = v_e)
      
      
      res_out = vcfac*(1 - res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]) +
        (1-vcfac)*(1 - res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]/
                     (res$i_pos_v_neg_tau3[seroout$study_years[[index_survey]]-1939,] +
                        res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]))
      
      seroprev_predict_survey_r = rbind(seroprev_predict_survey_r,  res_out)
      
      
    } # end of for(index_survey)
    
    seroprev_predict_survey_print_r = data.frame(sero_studies=rep(seroout$sero_studies, times = 1),
                                                 seroprev_predict_survey_r)
    
    
    if (i==1) {
      r_sero_predictions_lo=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
    }
    if (i==2) {
      r_sero_predictions=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
    }
    if (i==3) {
      r_sero_predictions_hi=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
    }
    
    
  } #end of for(i)
  
  
  ### PLOT ###
  par(mfrow=c(6,7),  oma = c(4, 4, 0, 0), mar = c(2, 2, 1, 1))
  for (i in 1:seroout$no_sero_surveys){
    #get upper bound for plot axes
    maxsero = max(as.numeric(r_sero_predictions_hi[i,1:86]),
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
    polygon(c(85:0,0:85), c(rev(as.numeric(r_sero_predictions_hi[i,1:86])),
                            (as.numeric(r_sero_predictions_lo[i,1:86]))) ,
            border=NA, col = rgb(1, 0, 0,0.2) )
    
    #add the median predictions
    lines(0:85,as.numeric(r_sero_predictions[i,1:86]),
          type="l",col = "red")
    
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
                   glm_mcmc_out, 
                   model_type = "Foi"){
  
  if(!anyNA(mcmc_out)){
    ### hpd for serology ###
    #full uncertainty for shared parameters
    hpd_sero = hpd_out(mcmc_out[,grep("^vac", names(mcmc_out))])
    
    if(model_type == "Foi"){
      #hpd by model type
      hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out[,grep("^Foi", names(mcmc_out))]) )
    } else {
      hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out[,grep("^R0", names(mcmc_out))]) )
    }
    
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
  
  
  if(model_type == "Foi"){
    df <-  data.frame(adm0 = shp1$GID_0,
                      adm0_adm1 = shp1$GID_1,
                      FOI = shp1$runs)
    
    image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
               axis.args = list(at = c(-7:1), labels =c("1e-07", "1e-06","1e-05","1e-04", "0.001", "0.01", "0.1","1","10"), las =2),
               legend.mar = 3.5)
  } else  {
    df <-  data.frame(adm0 = shp1$GID_0,
                      adm0_adm1 = shp1$GID_1,
                      R0 = shp1$runs+1)
    image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
               axis.args = list(at = c(-3:1), labels =c("1.001", "1.01", "1.1","2","11"), las =2),
               legend.mar = 3.5)
  }
  
  return(df)
  
}

#------------------------------------------------------------------------------------------------------#
plot_prior_post2 = function(mcmc_out,
                            prior_col = "red",
                            glm_or_sero ){
  
  if(glm_or_sero == "GLM"){
    #set up output plot
    par(mfrow=c(3,2), mar = c(4,4,2,1) + 0.1)
    
    #get the variables that are country factors
    jj = grep("adm05", names(mcmc_out))
    sd.prior = 2 #hardcoded- check matches GLMprior
    
    p=as.data.frame(mcmc_out)[,names(mcmc_out)[jj]]
    t = as.vector(as.matrix(p) )
    
    #plot for country factors
    plot(density(t),
         ylim=c(0,1), main="", xlab = "Country factor")
    vec = seq(min(t)-100,max(t)+100, by = (max(t)-min(t))/1000)
    polygon(vec, dnorm(vec, mean = 0, sd = sd.prior), col = prior_col, border = prior_col)
    polygon(density(t), col = rgb(0,0,0,0.5))
    
    
    #plot for everything else
    kk = grep("adm05", names(mcmc_out), invert = TRUE)[1:10]
    for (i in kk){
      q = as.numeric(as.data.frame(mcmc_out)[,names(mcmc_out)[i]])
      
      plot(density(q), main = "", xlab=names(mcmc_out)[i])
      vec = seq(min(q)-100, max(q)+100, by = (max(q)-min(q))/1000)
      polygon(vec,dnorm(vec, mean = 0, sd = 30), col = prior_col, border = prior_col)
      polygon(density(q), col = rgb(0,0,0,0.5))
    }
    
  }
  
  if(glm_or_sero == "SERO"){
    
    par(mfrow=c(1,2))
    
    #plot vaccine efficacy
    plot(density(exp(mcmc_out$vac_eff)), xlim = c(min(exp(mcmc_out$vac_eff)),1),
         main = "", xlab="Vaccine efficacy")
    
    vec = seq(min(exp(mcmc_out$vac_eff)), 1.2, length.out = 1000)
    polygon(vec, dnorm(vec, mean = 0.975, sd = 0.05), col = prior_col, border = prior_col)
    
    polygon(density(exp(mcmc_out$vac_eff)), col = rgb(0,0,0,0.5))
    
    
    #plot vaccine coverage factor for CMRs
    plot(density(exp(mcmc_out$vc_factor_CMRs)), xlim = c(0,1),
         main = "", xlab="Vaccine factor CMRs")
    
    vec = seq(-0.01, 1.01, length.out = 1000)
    polygon(vec,dunif(vec), col = prior_col, border = prior_col)
    
    polygon(density(exp(mcmc_out$vc_factor_CMRs)), col = rgb(0,0,0,0.5))
  }
}

#------------------------------------------------------------------------------------------------------#

foi_prevac2 = function(adm,
                       R0,
                       pop_moments,
                       polydeg = 6) {
  
  lambda = NULL
  if (polydeg > 0) {
    for (deg in 1:polydeg) {
      if (is.na(max(adm)) == 0) {
        PM = as.numeric(pop_moments[, deg + 1])
      } else {
        PM = pop_moments[deg]
      }  #INDEX IS DIFFERENT AS INCLUDES NAMES
      lambda = cbind(lambda, (-1)^(deg + 1) * PM/factorial(deg - 1))  #
    }
  }
  lambda[, 1] = lambda[, 1] - 1/R0  # we have to put the equation =0, so the term of order 0 (first column) should integrate -1/R0
  
  out = sapply(1:nrow(lambda), function(j) polyroot(lambda[j, ]))
  
  out[abs(Arg(out)) <= 1e-10] = Re(out)[abs(Arg(out)) <= 1e-10]
  out[abs(Arg(out)) > 1e-10] = NA  # here we have a resolution problem : in case of survey=5 and R0=5 (for instance), we have no solution in Real
  dt = dim(out)
  out = as.numeric(out)
  dim(out) = dt
  if (polydeg > 2) {
    out = apply(out, 2, min, na.rm = T)
  }
  if (is.na(max(adm)) == 0){
    names(out) = paste("FOI", pop_moments[,1], sep = "_")
  }  #NAME IS ATTACHED TO POP_MOMENTS FOR WHOLE
  return(out)
}

#------------------------------------------------------------------------------------------------------#

create_R0_lookup2 = function(dim_adm,
                             dat,
                             dim_year,
                             dim_age,
                             p_prop_3d,
                             P_tot_2d,
                             inc_v3d,
                             pop_moments_whole,
                             pop_moments_agg,
                             vac_eff_arg = 1,
                             granularity = "coarse") {
  
  dn_R0 = switch(granularity,
                 "fine" = c(seq(1, 1.099 , by=0.001), seq(1.1, 1.198, b=0.002), seq(1.2, 1.495, by=0.005), seq(1.5, 1.99, by=0.01),
                            seq(2, 2.9, by= 0.05 ), seq(3,3.9, by = 0.1), seq(4,10, by=0.2)),
                 "medium" = c(seq(1,2,by=0.01), seq(2.1,4.9,by=0.1), seq(5,10, by=1) ),
                 "coarse" = c(seq(1,1.9,by=0.1), seq(2,4.5,by=0.5), seq(5,10, by=1) ) )
  
  
  R0_lookup = rep(NA, length(dim_adm)*length(dn_R0))
  dim(R0_lookup) = c(length(dim_adm),length(dn_R0))
  colnames(R0_lookup) = dn_R0
  rownames(R0_lookup) = dat$adm0_adm1
  adm_ind = 1:length(dim_adm)
  
  ptm = proc.time()
  for (r in 1:length(dn_R0)){
    print(paste(r,"/", length(dn_R0), ";   R0=", dn_R0[r]))
    R0_rep = rep(dn_R0[r], length(dim_adm))
    
    tmp = R0_recurrence_seroprev_whole2(adm = adm_ind,
                                        dat = dat,
                                        R0 = R0_rep,
                                        t0_vac_adm = t0_vac_africa,
                                        dim_year = dim_year,
                                        dim_age = dim_age,
                                        p_prop_3d = p_prop_3d,
                                        P_tot_2d = P_tot_2d,
                                        inc_v3d = inc_v3d,
                                        pop_moments_whole = pop_moments_whole,
                                        vac_eff_arg = vac_eff_arg
    )$Ninf_t_province # this is for each year from 1940 to 2050
    
    R0_lookup[,r] = rowSums(tmp[, which(dim_year==1984):which(dim_year==2019)])
  }
  save(R0_lookup, file="R0_lookup_table.Rdata")
  
}

#------------------------------------------------------------------------------------------------------#
R0_recurrence_seroprev_whole2 = function(adm,
                                         dat,
                                         R0,
                                         t0_vac_adm,
                                         dim_year,
                                         dim_age,
                                         p_prop_3d,
                                         P_tot_2d,
                                         inc_v3d,
                                         pop_moments_whole,
                                         vac_eff_arg) {
  #
  #THIS NOW ACTS AS A WRAPPER FOR fun_recurrence_seroprev WHICH CONTAINS THE one COPY OF THE EQUATIONS
  #COMBINING THE FUNCTIONS IN THIS WAY SHOULD REDUCE INCONSISTENCIES BETWEEN DIFFERENT MODEL VERSIONs
  #SUCH AS: THE INCLUSION OF ENVIRONMENTAL TRANSMISSION
  
  # Declare outputs
  n_inf_at = s_atau3 = s_atau2= s_atau1= Stot_at = Ninf_t = NULL
  
  # create the 4 compartments accounting for natural infection i (+/-) and vaccination v (+/-)
  i_neg_v_neg_tau3 = i_pos_v_neg_tau3 = i_neg_v_pos_tau3 = i_pos_v_pos_tau3 = NULL
  
  # I firstly need to calculate for each adm the FOI at t0_vac
  pop_mom=NULL
  for(j in 1:length(adm)){
    pop_mom = rbind (pop_mom, pop_moments_whole[adm[j], ])
  }
  pop_mom = cbind(dat$adm0_adm1[adm], pop_mom)
  
  #get foi before vaccination
  foi = YFburden::foi_prevac(adm = adm, R0=R0, pop_moments = pop_mom, polydeg=6)
  
  for (i in 1:length(adm) ){
    
    time_limit = min(1981, max(t0_vac_adm[i], 1940) )  
    p_at= p_prop_3d[adm[i], ,]
    P_tot=P_tot_2d[adm[i], ]
    inc_v=inc_v3d[adm[i], ,]
    
    seroprev_out = R0_recurrence_seroprev(R0 = R0[i],
                                          foi = foi[i],
                                          foi_const = 0,
                                          time_limit = time_limit,
                                          dim_year = dim_year,
                                          dim_age = dim_age,
                                          p_at = p_at,
                                          P_tot = P_tot,
                                          inc_v = inc_v,
                                          vac_eff_arg = vac_eff_arg
    )
    
    if(i==1){
      
      i_neg_v_neg_tau3 = seroprev_out$i_neg_v_neg_tau3 #FIRST JUST ASSIGN MATRIX
      i_pos_v_neg_tau3 = seroprev_out$i_pos_v_neg_tau3
      i_neg_v_pos_tau3 = seroprev_out$i_neg_v_pos_tau3
      i_pos_v_pos_tau3 = seroprev_out$i_pos_v_pos_tau3
      Ninf_t = seroprev_out$Ninf_t
      Stot_at = seroprev_out$Stot
      n_inf_at = seroprev_out$n_inf
      
    } else if(i==2){
      
      i_neg_v_neg_tau3 = abind::abind(i_neg_v_neg_tau3, seroprev_out$i_neg_v_neg_tau3, along=0) #THEN BIND ALONG NEW FIRST DIMENSION, RESULT IS AN ARRAY
      i_pos_v_neg_tau3 = abind::abind(i_pos_v_neg_tau3, seroprev_out$i_pos_v_neg_tau3, along=0)
      i_neg_v_pos_tau3 = abind::abind(i_neg_v_pos_tau3, seroprev_out$i_neg_v_pos_tau3, along=0)
      i_pos_v_pos_tau3 = abind::abind(i_pos_v_pos_tau3, seroprev_out$i_pos_v_pos_tau3, along=0)
      Ninf_t = abind::abind(Ninf_t, seroprev_out$Ninf_t,along=0)
      Stot_at = abind::abind(Stot_at, seroprev_out$Stot,along=0)
      n_inf_at = abind::abind(n_inf_at, seroprev_out$n_inf,along=0)
      
    } else {
      
      i_neg_v_neg_tau3 = abind::abind(i_neg_v_neg_tau3, seroprev_out$i_neg_v_neg_tau3, along=1) #NOW FIRST DIMENSION EXISTS, BIND ALONG THAT
      i_pos_v_neg_tau3 = abind::abind(i_pos_v_neg_tau3, seroprev_out$i_pos_v_neg_tau3, along=1)
      i_neg_v_pos_tau3 = abind::abind(i_neg_v_pos_tau3, seroprev_out$i_neg_v_pos_tau3, along=1)
      i_pos_v_pos_tau3 = abind::abind(i_pos_v_pos_tau3, seroprev_out$i_pos_v_pos_tau3, along=1)
      Ninf_t = abind::abind(Ninf_t, seroprev_out$Ninf_t,along=1)
      Stot_at = abind::abind(Stot_at, seroprev_out$Stot,along=1)
      n_inf_at = abind::abind(n_inf_at, seroprev_out$n_inf,along=1)
      
    }
    
  }
  
  return(list(i_neg_v_neg_tau3 = i_neg_v_neg_tau3,
              i_pos_v_neg_tau3 = i_pos_v_neg_tau3,
              i_neg_v_pos_tau3 = i_neg_v_pos_tau3,
              i_pos_v_pos_tau3 = i_pos_v_pos_tau3,
              Stot_at = Stot_at,
              Ninf_t_province = Ninf_t,
              n_inf_at_province = n_inf_at))
}

#------------------------------------------------------------------------------------------------------#

fun_calc_transmission_Africa2 = function(x,
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
                                  polydeg=5,
                                  R0_lookup,
                                  model_type) {
  
  #get aggregated vc and pop over observation period
  aggout=create_pop30_agg_vc30_agg(pop1, vc2d)
  
  #glm predictions
  mypreds_nc  = fun_calcPred(coefs = as.numeric(params)[ii],
                             newdata=x,
                             type="link",
                             varsin=varsin_nc)
  
  #probability of detection
  p_detect =  fun_calc_pdetect_multi_both2(x,
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
                                           pop_moments_whole,
                                           varsin_nc,
                                           aggout$vc30_agg,
                                           aggout$pop30_agg,
                                           model_type)
  
  
  p_detect_link = mean(p_detect)
  
  #calculating number of infections over the observation period for the whole region
  Ninf_whole = exp( mypreds_nc - p_detect_link)
  
  #translate this number into a transmission intensity for each model type
  if(model_type == "Foi"){
    pop_vc_moments = aggout$pop_vc_moments
    
    if(polydeg>ncol(pop_vc_moments)) error("fun_calc_transmission_Africa: invalid value for polydeg.\n")
    
    z = -Ninf_whole
    
    if(polydeg>0) for(i in 1:polydeg) {
      z = cbind(z,(-1)^(i+1)*pop_vc_moments[,i+1]/factorial(i-1))
    }
    
    transmission_whole = sapply(1:nrow(x), function(i) polyroot(z[i,]))
    transmission_whole[abs(Arg(transmission_whole))<=1e-10] = Re(transmission_whole)[abs(Arg(transmission_whole))<=1e-10]
    transmission_whole[abs(Arg(transmission_whole))>1e-10] = NA
    
    dt = dim(transmission_whole)
    transmission_whole = as.numeric(transmission_whole)
    dim(transmission_whole) = dt
    transmission_whole = apply(transmission_whole,2,min,na.rm=T)
    
  } else if(model_type == "R0"){
    
    transmission_whole = rep(NA, length(Ninf_whole))
    
    for (i in 1: length(t0_vac_africa) ){
      inf_bound = findInterval(Ninf_whole[i], R0_lookup[i,])
      sup_bound = inf_bound + 1
      
      if( sup_bound <= ncol(R0_lookup)){
        x_R0_1 = R0_lookup[i,inf_bound]# find inf bound
        x_R0_2 = R0_lookup[i,sup_bound] # sup bound
        
        y_R0_1 = as.numeric(colnames(R0_lookup)[inf_bound]) # find inf bound
        y_R0_2 = as.numeric(colnames(R0_lookup)[sup_bound]) # sup bound
        
        ## manually calculate the linear interpolation
        a_lin_inter = (y_R0_2 - y_R0_1)/(x_R0_2-x_R0_1)
        b_lin_inter = y_R0_1 - a_lin_inter*x_R0_1
        
        transmission_whole[i] = a_lin_inter*Ninf_whole[i] + b_lin_inter
        
      } else if (sup_bound  > ncol(R0_lookup)){
        transmission_whole[i] = as.numeric(colnames(R0_lookup)[ncol(R0_lookup)])
      }
    }
    
  }
  
  return(transmission_whole)
}


#------------------------------------------------------------------------------------------------------#
fun_calc_pdetect_multi_both2 = function(x,
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
                                        pop_moments_whole,
                                        varsin_nc,
                                        vc30_agg,
                                        pop30_agg,
                                        model_type){
  
  adm1s = seroout$adm1s
  sero_studies = seroout$sero_studies
  no_sero_surveys = seroout$no_sero_surveys
  
  # calculate the b term from the left part of equation 6(plos Med paper)
  adm = match(unlist(adm1s), vc30_agg$adm0_adm1)
  t0_vac_adm = t0_vac_africa[adm]
  
  # now I have to calculate the Number of infections over observation period among each PROVINCE covered by a survey
  vac_eff = params[1]
  
  if(model_type == "R0"){
    
    R0_vec = NULL
    for(k in 1:no_sero_surveys) { # repeat the numbers of time of each adm1s
      R0_vec = c(R0_vec,
                 rep(params[which(names(params)==paste0("R0_", sero_studies[k]))],
                     length(adm1s[[k]])) )
    }
    
    res = R0_recurrence_seroprev_whole(adm = adm,
                                       dat = dat,
                                       R0 = R0_vec,
                                       t0_vac_adm,
                                       dim_year,
                                       dim_age,
                                       p_prop_3d,
                                       P_tot_2d,
                                       inc_v3d,
                                       pop_moments_whole,
                                       vac_eff_arg = vac_eff)
    
    # res returns the total Nb of infection on the whole period, I only need 1984 to 2019
    Ninf = rowSums( res$Ninf_t_province[,which(dim_year==1984):which(dim_year==2019)] )
    Ninf = ifelse(Ninf<1, 1, Ninf)
    
  } else if(model_type == "Foi"){
    
    vec_foi = vec_vc_factor = NULL
    for(i in 1:no_sero_surveys) {
      vec_foi = c(vec_foi,
                  rep(params[which(names(params)==paste("Foi_",sero_studies[i],sep=""))],
                      length(adm1s[[i]])))
      if(!is.na(seroout$vc_factor[i])) {
        vec_vc_factor = c(vec_vc_factor, rep(seroout$vc_factor[i],length(adm1s[[i]])))
      } else {
        vec_vc_factor = c(vec_vc_factor,
                          rep(params[which(names(params)==paste("vc_factor_",sero_studies[i],sep=""))],
                              length(adm1s[[i]])))
      }
    }
    
    vec_foi = vec_foi[match(vc30_agg[adm,1],unlist(adm1s))]
    vec_vc_factor = vec_vc_factor[match(vc30_agg[adm,1], unlist(adm1s))]
    
    popvc = (1-vac_eff*vec_vc_factor*vc30_agg[adm,-(1)])*pop30_agg[adm,-(1)]
    expfoi = exp(outer(-vec_foi, 0:100))
    Ninf = vec_foi*rowSums(expfoi*popvc, na.rm = T)
  }
  
  #calcualte glm predictions
  mypreds_nc = fun_calcPred(coefs = params[ii],
                            newdata=x,
                            type="link",
                            varsin=varsin_nc)
  
  #calculate b
  b = mypreds_nc[adm]-log(Ninf)
  names(b) = unlist(adm1s)
  
  return( b )
}
