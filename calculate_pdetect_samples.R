
# script to calcualte p_detect for Asia

library(mcmcplots)
library(ggmcmc)
library(gridExtra)
library(maptools)
library(readr)
library(fields)
library(Hmisc)
library(viridis)
library(R.utils)
library(magrittr)
library(ggplot2)
library(dplyr)

library(YFestimation)
library(snapalette)
library(KsetupR)

sourceDirectory("Functions", modifiedOnly = FALSE)

filepath_to_use <- readRDS("sero_filepath_to_use.RDS")

if(!exists("mcmc_out")){
  mcmc_out <- get_chains(paste0("../YellowFeverModelEstimation2019/", filepath_to_use),
                         burnin = 1, thin = 1)
  
  write.csv(sample(exp(mcmc_out$vac_eff), 1000), "vac_eff_samples.csv", row.names = FALSE)
}

glm_mcmc_out <- YFestimation::get_chains(paste0("GLM_MCMC_chain_20200124_bestglm_6_", 
                                                1),
                                         burnin =6e4, thin = 100)

names(glm_mcmc_out) = gsub("^log.", "", names(glm_mcmc_out))

### LOADING SEROLOGY DATA ###
Serology = read.csv(paste0("../Data/","Serology/Serology_newgadm_2019.csv"), 
                    stringsAsFactors = FALSE)

Serology$gadm36 = gsub("_", "." ,Serology$gadm36)

seroout = process_serology(Serology, adm = "gadm36")
for(s in 1:seroout$no_sero_surveys){
  seroout$adm1s[[s]] = paste(seroout$adm1s[[s]], "_1", sep = "")
}

for(model_ind in 1:20){
  
  #remove families of NHP that are not to be included
  model_form_whole = read.csv("bestglm_6.csv", stringsAsFactors = FALSE) %>% dplyr::select(-Criterion)
  
  ### POPULATION AND VACCINATION DATA ###
  if(!exists("all_res_pop_3d")){
    all_res_pop_3d = get_pop_data_3d() 
  }
  pop1 = all_res_pop_3d$pop1                                      
  P_tot = all_res_pop_3d$P_tot                                   
  p_prop = all_res_pop_3d$p_prop                                   
  
  pop1 %<>% rename(adm0_adm1 = adm1)
  
  pop1 %<>% arrange(year) %>% dplyr::select(-c(country_code, country))
  
  names(pop1)[3:ncol(pop1)] = paste0("a", 0:100)
  
  pop1 %<>% filter(year<2051)
  p_prop %<>% filter(year<2051)
  P_tot %<>% filter(year<2051)
  
  # VACCINATION DATA #
  vc2d = readRDS("../YellowFeverModelEstimation2019/vacc_1940_1950.RDS")
  vc2d = as.data.frame(vc2d)
  
  # expand vc2d to match pop1
  vc2d %<>% mutate(year = as.numeric(as.character(year)))
  vc_tmp = pop1 %>% filter(!adm0_adm1 %in% vc2d$adm0_adm1)
  vc_tmp[, grep("adm0_adm1|year", names(vc_tmp), invert = TRUE)] = 0
  
  vc2d %<>% bind_rows(vc_tmp)
  vc2d %<>% unique()
  vc2d %<>% filter(!is.na(adm0_adm1))
  
  
  
  
  covar_family = names(model_form_whole)[ which(model_form_whole[model_ind, ] == TRUE) ] 
  
  # LOAD ENVIRONMENTAL DATA #
  Env_Table_path <- "../Data/Environment/global_dat"
  
  filename <- KsetupR::get_latest_file(path = Env_Table_path, pattern = "dat_wes_mosquito")
  
  dat <- read.csv(filename, stringsAsFactors = FALSE)
  
  dat <- dat[, -which(!names(dat) %in% covar_family & grepl("family", names(dat)))] 
  
  
  glm_mcmc_out <- YFestimation::get_chains(paste0("GLM_MCMC_chain_20200124_bestglm_6_", 
                                                  model_ind),
                                           burnin =6e4, thin = 100)
  
  names(glm_mcmc_out) = gsub("^log.", "", names(glm_mcmc_out))
  
  # get model covariates
  covar = names(glm_mcmc_out %>% dplyr::select(-c(Intercept, 
                                                  starts_with("adm05"), 
                                                  posteriorProb, 
                                                  acceptRate))) 
  
  
  
  # make extra elements and normalise
  dat %<>% dplyr::select(-year)
  dat %<>% filter( is.finite(MIR.max))
  dat %<>% filter(!is.na(temp_mean))
  #dat %<>% tidyr::drop_na()
  dat %<>% mutate(continent = ifelse(adm0 %in% c("CIV", "STP"),
                                     "Africa", continent))
  dat = adjust_env_dat(dat)
  
  dat %<>% mutate(adm0_adm1 = adm1)
  
  
  
  
  #-------------------------------------------------------
  # filter vac and pop to match dat
  vc2d %<>% filter(adm0_adm1 %in% dat$adm1)
  pop1 %<>% filter(adm0_adm1 %in% dat$adm1)
  
  
  #-------------------------------------------------------
  #vc3d
  vc3d = transform_into_vc3d(vc2d) 
  #------------------------------------------------------
  class(vc3d) = "numeric"
  # t0_vac_africa #
  t0_vac_africa = calc_t0_vac_africa(vc3d)
  
  # inc_v3d #
  if(!exists("inc_v3d")){
    inc_v3d = calc_incidence_vac_general(vc3d)
  }
  
  # CALCULATE population moments #
  p_prop_3d = array( data = p_prop$population,
                     dim = c(length(unique(p_prop$adm1)),
                             length(unique(p_prop$year)),
                             length(unique(p_prop$age))),
                     dimnames = list(unique(p_prop$adm1),
                                     unique(p_prop$year),
                                     unique(p_prop$age)))
  
  pop_moments_whole = calc_pop_moments(p_prop_3d, 
                                       t0_vac_africa,
                                       dim_adm = dimnames(vc3d)[[1]],
                                       dim_year = dimnames(vc3d)[[2]],
                                       dim_age = dimnames(vc3d)[[3]])
  
  #aggregate
  if(!exists("list_aggregate_pop_vc")){
    list_aggregate_pop_vc = 
      YFestimation::Make_aggregate_pop_vc_3d(pop1=pop1, 
                                             vc2d=vc2d %>% dplyr::select(-adm0), 
                                             sero_studies=seroout$sero_studies, 
                                             adm1s=seroout$adm1s) 
  }
  pop_agg3d = list_aggregate_pop_vc$pop_agg3d 
  vc_agg3d = list_aggregate_pop_vc$vc_agg3d 
  
  pop_agg3d[is.na(pop_agg3d)] = 0
  vc_agg3d[is.na(vc_agg3d)] = 0
  
  #calculate aggregated incidence (same function as before)
  class(vc_agg3d) = "numeric"
  inc_v3d_agg = calc_incidence_vac_general(vc_agg3d)  
  
  #calculate aggregated moments (different fucntion before)
  class(pop_agg3d) = "numeric"
  pop_moments_agg = calc_pop_moments_agg(pop_agg3d,
                                         seroout$t0_vac,
                                         dim_year = unique(pop1$year),
                                         seroout$study_years) 
  
  
  
  # CREATE R0 LOOKUP TABLE #
  load(paste0("../YellowFeverModelEstimation2017/",
              "R0_lookup_table.Rdata") )
  
  
  # pop at survey #
  foi_const_surv = rep(0, seroout$no_sero_surveys)
  
  list_pop_at_survey = create_pop_at_survey(pop_agg3d, 
                                            seroout$sero_studies, 
                                            dim_year = unique(pop1$year))
  p_at_survey = list_pop_at_survey$p_at_survey_3d
  P_tot_survey = list_pop_at_survey$P_tot_survey_2d
  
  model_form <- paste0("cases_or_outbreaks~", paste(covar, collapse = "+")  ,  "+adm05") 
  
  
  object_glm <- YFestimation::fit_glm(dat = dat, 
                                      depi = match("cases_or_outbreaks", 
                                                   names(dat)), 
                                      models = model_form)
  
  beta0 <- object_glm[[1]]
  x <- object_glm[[2]]
  y <- object_glm[[3]]
  
  
  
  
  
  out = get_hpd(mcmc_out, glm_mcmc_out)
  
  Foi_param = out$Foi_param
  
  ii = c(2:(grep("Foi", names(Foi_param[,1]))[1]-1))
  varsin_nc = grep("adm05|surv|posterior|acc|dtp", names(glm_mcmc_out), invert = TRUE)
  
  colours = snapalette("Pop", 100, type = "continuous")
  
  vc2d %<>% arrange( adm0_adm1)
  
  #order x by FOI_param
  x = x[, names(Foi_param[ii,2])]
  
  n_samples = 10
  pdetect_out = NULL
  
  for(i in 1:n_samples){
    sample_ind1 = sample(1:nrow(mcmc_out),1)
    sample_ind2 = sample(1:nrow(glm_mcmc_out),1)
    
    parm = bind_cols(vac_eff = exp(mcmc_out[sample_ind1, "vac_eff"]),
                     glm_mcmc_out[sample_ind2, 1:(ncol(glm_mcmc_out)-2)],
                     exp(mcmc_out[sample_ind1, grep("^Foi", names(mcmc_out))]),
                     vc_factor_CMRs = exp(mcmc_out[sample_ind1, "vc_factor_CMRs"]))
    class(parm) = "numeric"
    
    
    #get aggregated vc and pop over observation period
    aggout = create_pop30_agg_vc30_agg(pop1, vc2d %>% select(-adm0))
    
    #glm predictions
    mypreds_nc  = fun_calcPred(coefs = as.numeric(parm)[ii],
                               newdata=x,
                               type="link",
                               varsin=varsin_nc)
    
    #probability of detection
    p_detect =  fun_calc_pdetect_multi_both(x,
                                            ii,
                                            seroout,
                                            parm,
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
                                            model_type = "Foi")
    
    
    p_detect_link = mean(p_detect)
    
    pdetect_out %<>% rbind(p_detect_link)
  }
  
  saveRDS(pdetect_out, paste0("6" ,"_pdetect_", model_ind, ".RDS"))
}
