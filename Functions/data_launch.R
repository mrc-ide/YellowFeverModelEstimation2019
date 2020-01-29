# Launch data #

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

adjust_env_dat = function(dat) {
  
  #which presence/abence outcome to consider
  depvar = "cases_or_outbreaks"
  depi = match(depvar, names(dat)) # column number for the chosen outcome
  
  
  # adding a categorical "dominant land cover" variable:
  LC_i = grep("LC",names(dat))
  LC_dom = names(dat)[LC_i][  as.numeric(apply(dat[,LC_i],1, which.max ))]
  dat = cbind(dat, LC_dom) # 37 potential covariates/levels
  dat$LC_dom = as.numeric(as.factor(dat$LC_dom))
  
  #setting NA to 0
  dat$surv.qual.adm0[is.na(dat$surv.qual.adm0)] = 0
  
  # setting the surv.qual.adm0 for AGO to 0 (very few cases reported):
  dat$surv.qual.adm0[dat$adm0=="AGO"] = 0
  
  
  #adding log(surv_qual_adm0)
  dat %<>% mutate(log.surv.qual.adm0 = ifelse(is.infinite(log(surv.qual.adm0)),
                                              0, log(surv.qual.adm0)))
  
  # a same categorical variable for all countries within the YFSD
  # aslo group countries not considered high or moderate risk together
  
  dat = dplyr::mutate(dat,
                      adm05 = ifelse(surv.qual.adm0>0,
                                     "AFR",
                                     ifelse(risk %in% c("high", "moderate"),
                                            adm0,
                                            "low_risk")))
                                     # 
  #add continental factor
  dat = dplyr::mutate(dat,
                      continent_factor = ifelse(continent=="Africa",
                                                1,
                                                0))
  
  # aggregate monkeys
  dat = dplyr::mutate(dat,
                      aggregate_family = rowSums(dat[, grep("family", names(dat))]))
  
  # aggregate monkeys
  dat = dplyr::mutate(dat,
                      aggregate_genus = rowSums(dat[, grep("genus", names(dat))]))
  
  
  
  # add ranges
  dat = dplyr::mutate(dat,
                      EVI.range = EVI.max - EVI.min,
                      temp_range = temp_max - temp_min)
  
  
  # add simple temp suitability
  parm = read.csv("temp_suit_parm.csv", stringsAsFactors = F)
  parm2 = as.numeric(parm$x)
  names(parm2) = parm$X
  dat = dplyr::mutate(dat,
                      temp_suit_mean = temp_suitability(dat$temp_mean, parm2),
                      temp_suit_min = temp_suitability(dat$temp_min, parm2),
                      temp_suit_max = temp_suitability(dat$temp_max, parm2))
  
  # sorting risk
  dat = dplyr::mutate(dat, risk = as.numeric(as.factor(risk)))
  dat$risk[is.na(dat$risk)] = max(dat$risk, na.rm = TRUE)+1
  
  # add DTP
  dtp = read.csv("../Data/Healthcare/Outputs/dtp3_coverage_by_country.csv", stringsAsFactors = FALSE)
  dtp %<>% rename(adm0 = iso)
  dat %<>% left_join(dtp %>% dplyr::select(-name))
  dat$dtp3_coverage[is.na(dat$dtp3_coverage)] = 0
  
  v1 = apply(dat,2,var, na.rm = TRUE)
  for(i in 9:(ncol(dat))) {
    
    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {
      
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
  
  
  pop1 = read.csv("Z:/Data/Population/all_pop_at_adm1_1940-2100_landscan2017_gadm36.csv",
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
# ------------------------------------------------------------------------------------------------------------------------------------------------------------

fit_glm2 = function(dat, depi, models) {
  
  fm_best = models[1] #take first model
  
  bm = glm(as.formula(fm_best), data=dat, family=binomial(link="cloglog"))
  
  # setting up the evaluation of the likelihood:
  beta = coefficients(bm)
  
  vl=NULL
  for(i in 1:ncol(dat)) {
    if(length(grep(names(dat)[i], names(beta)))>0) {
      vl = c(vl,i) # select the variables used in the GLM
    }
  }
  
  
  x = cbind(Intercept=1,dat[,vl])
  j_expand = !sapply(1:ncol(x), function(i) is.numeric(x[,i]) & !is.factor(x[,i]))
  
  # create indicative variables for adm05
  x %<>% mutate(adm05 = as.factor(paste0("adm05", adm05)))
  x2 = tidyr::spread(x, adm05, adm05)
  for(i in grep("adm05", names(x2))){
    x2[, i] = as.numeric(x2[, i])
  }
  x2[is.na(x2)] = 0
  
  y = dat[,depi]  # dependant variable vector
  
  beta0 = coefficients(bm)
  names(beta0)[names(beta0)=="(Intercept)"] = "Intercept"
  
  mm = match(names(beta0),colnames(x2))
  x2 = x2[,mm]
  
  beta0 = beta0[match(colnames(x2),names(beta0))]
  
  return(list(beta0=beta0, x=x2, y=y))
  
}
