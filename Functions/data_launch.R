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
  # adding a categorical "dominant land cover" variable:
  LC_i = grep("LC",names(dat))
  LC_dom = names(dat)[LC_i][ apply(dat[,LC_i],1, which.max )]
  dat = cbind(dat, LC_dom) # 37 potential covariates/levels
  
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
  
  #add continental factor
  dat = dplyr::mutate(dat,
                      continent_factor = ifelse(continent=="Africa",
                                                1,
                                                0))

  # If we fit the model on dat and if we want to project estimates on dat, variable from dat
  # need to be expressed on the same scale than those from dat , thus we normalize dat relatively to dat
  v1 = apply(dat,2,var, na.rm = TRUE)
  for(i in 9:(ncol(dat))) {
    
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


# ------------------------------------------------------------------------------------------------------------------------------------------------------------
get_pop_data_3d = function() {
  
  
  pop1 = read.csv("../Data/Population/pop_at_adm1_1940-2100_landscan2017_gadm36.csv",
                  stringsAsFactors = FALSE)

  pop2d = tidyr::gather(pop1, age, population, -c(country_code, adm1, country, year))
  
  out = Make_Ptot_Pprop(pop2d)
  
  return(list(pop1 = pop1,  P_tot = out$P_tot, p_prop = out$p_prop))
}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
Make_Ptot_Pprop = function(pop2d){
  
  # get totals
  P_tot = pop2d %>% dplyr::group_by(adm1, year) %>% 
    dplyr::summarise(pop = sum(population, na.rm = TRUE)) %>%
    unique()
  
  #get proportions
  p_prop = pop2d %>% left_join(P_tot, by = c("adm1", "year"))
  p_prop %<>% 
    dplyr::mutate(population = population/pop) %>% 
    dplyr::select(-pop)
  
  return(list( P_tot = P_tot, p_prop = p_prop))
  
}
