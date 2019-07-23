# Launch data #

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

adjust_env_dat = function(dat) {
  
  #which presence/abence outcome to consider
  depvar = "cases_or_outbreaks"
  depi = match(depvar, names(dat)) # column number for the chosen outcome
  
  # # excluding LC1,3,5,15 as they only ever cover <5% -
  # ex_i = match(c("LC1","LC3","LC5","LC15"), names(dat))
  # dat = dat[,-ex_i]
  # 
  # # adding a categorical "dominant land cover" variable:
  # LC_i = grep("LC",names(dat))
  # LC_dom = names(dat)[LC_i][ apply(dat[,LC_i],1, which.max )]
  # dat = cbind(dat, LC_dom) # 37 potential covariates/levels
  
  #setting NA to 0
  dat$surv.qual.adm0[is.na(dat$surv.qual.adm0)] = 0
  
  # setting the surv.qual.adm0 for AGO to 0 (very few cases reported):
  dat$surv.qual.adm0[dat$adm0=="AGO"] = 0
  
  
  #adding log(surv_qual_adm0)
  dat %<>% mutate(log.surv.qual.adm0 = ifelse(is.infinite(log(surv.qual.adm0)),
                                              0, log(surv.qual.adm0)))
  
  # a same categorical variable for all countries within the YFSD
  dat = dplyr::mutate(dat,
                      adm05 = ifelse(surv.qual.adm0>0,
                                     "AFR", adm0))
  

  # If we fit the model on dat and if we want to project estimates on dat, variable from dat
  # need to be expressed on the same scale than those from dat , thus we normalize dat relatively to dat
  v1 = apply(dat,2,var, na.rm = TRUE)
  for(i in 8:(ncol(dat))) {
    
    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {
      
      dat[,i] = dat[,i]/sqrt(v1[i])
      dat[,i] = dat[,i]/sqrt(v1[i])
    }
  }
  
  dat = dat[,names(dat)!="surv_qual_adm0"]
  depi = match(depvar,names(dat))
  
  # I do that because pop is ordered that way and we need to match both pop and dat
  dat %<>% arrange(adm1)
  
  return(dat)
}