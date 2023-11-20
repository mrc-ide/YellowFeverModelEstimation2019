#transdimensional mcmc for estimating both r0 and foi models 


run_estimation = function(run_id=1){
  
  
  library(dplyr)
  library(readr)
  library(reshape2)
  library(abind)
  library(mvtnorm)
  library(truncdist)
  library(R.utils)
  library(magrittr)
  library(tidyr)
  library(Matrix)
  
  #homemade
  library(YFburden)
  library(YFestimation)
  library(KsetupR)
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # NEW FUNCTIONS #
  R.utils::sourceDirectory("Functions", modifiedOnly = FALSE)
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### LOADING SEROLOGY DATA ###
  
  Serology = read.csv(paste0("../Data/","Serology/Serology_newgadm_2019.csv"), stringsAsFactors = FALSE)
  
  Serology$gadm36 = gsub("_", "." ,Serology$gadm36)
  
  seroout = process_serology(Serology, adm = "gadm36")
  for(s in 1:seroout$no_sero_surveys){
    seroout$adm1s[[s]] = paste(seroout$adm1s[[s]], "_1", sep = "")
  }
  
  
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### POPULATION AND VACCINATION DATA ###
  
  all_res_pop_3d = get_pop_data_3d() 
  
  pop1 = all_res_pop_3d$pop1                                      
  P_tot = all_res_pop_3d$P_tot                                   #total populations for each adm and year
  p_prop = all_res_pop_3d$p_prop                                   #proportions of population
  
  pop1 %<>% rename(adm0_adm1 = adm1)
  
  pop1 %<>% arrange(year) %>% select(-c( country))
  
  names(pop1)[4:ncol(pop1)] = paste0("a", 0:100)
  
  pop1 %<>% filter(year<2051)
  
  pop1 %<>% rename(adm0 = country_code)
  
  #########################################################################################################
  ### VACCINATION DATA ###
  #########################################################################################################
  
  if(!file.exists("vacc_1940_1950.RDS")){
    vac_in = readRDS("../Data/Vaccination/Outputs/vac_coverage_1940_2050_all_gadm36.RDS")
    
    vc2d <-  vac_in %>% select(-c(doses, population)) %>% spread(age, coverage)
    
    names(vc2d)[3:ncol(vc2d)] = paste0("a", 0:100)
    
    vc2d %<>% tibble::add_column(adm0 = substr(vc2d$new_id, 1,3), .after = "year")                         #rename countries as adm0
    
    vc2d %<>% rename(adm0_adm1 = new_id)                        #renames adm1 as adm0_adm1
    
    # formally "repair_vc_data" from FOI model in Kevin's folder
    for (colIndex in 4:ncol(vc2d)){                                      #before 1995, we have NA values for those aged >75
      vc2d[,colIndex] = ifelse(is.na(vc2d[,colIndex]), vc2d[,colIndex-1], vc2d[,colIndex])
    }
    
    # restrict to lines in dat
    #vc2d %<>% filter(adm0_adm1 %in% dat$adm0_adm1)
    
    vc2d %<>% filter(!is.na(adm0))
    
    saveRDS(vc2d, "vacc_1940_1950.RDS")
  } else {
    vc2d = readRDS("vacc_1940_1950.RDS")
  }
  
  vc2d = as.data.frame(vc2d)
  vc2d %<>% filter(!is.na(adm0_adm1))
  
  vc2d %<>% filter(adm0_adm1 %in% pop1$adm0_adm1)
  vc2d %<>% ungroup()  %>% mutate(year = as.numeric(as.character(year)))
  
  #########################################################################################################
  ### AGGREGATE POPULATION AND VACCINATION DATA ###
  #########################################################################################################
  
  if(!file.exists("agg_pop_vc.RData")){ 
    
    #aggregate
    agg=Make_aggregate_pop_vc(select(pop1, -adm0), select(vc2d, -adm0), seroout$sero_studies, seroout$adm1s)
    
    pop_agg=agg$pop_agg
    vc_agg=agg$vc_agg
    
    dim_survey = seroout$sero_studies
    dim_year = as.numeric(unique(pop_agg$year))
    dim_age = names(pop1 %>% select(-c(adm0_adm1, year, adm0)))
    
    length_survey=length(dim_survey);   length_year=length(dim_year);   length_age=length(dim_age)
    
    ## pass in 3d
    pop_agg3d=rep(NA, nrow(pop_agg)*length_age)
    dim(pop_agg3d)=c(length_survey, length_year, length_age)
    vc_agg3d=rep(NA, nrow(vc_agg)*length_age)
    dim(vc_agg3d)=c(length_survey, length_year,length_age)
    
    for(ageIndex in 1:length_age) { # dim_age
      for(yearIndex in min(dim_year):max(dim_year)) {
        
        mm = match(pop_agg$adm0[pop_agg$year == yearIndex], dim_survey)
        pop_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex] = pop_agg[pop_agg$year==yearIndex, ageIndex+3]
        vc_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex] = vc_agg[vc_agg$year==yearIndex, ageIndex+3]
      }
    }
    
    dimnames(pop_agg3d)[[1]] = dimnames(vc_agg3d)[[1]] = dim_survey
    dimnames(pop_agg3d)[[2]] = dimnames(vc_agg3d)[[2]] = dim_year
    dimnames(pop_agg3d)[[3]] = dimnames(vc_agg3d)[[3]] = dim_age
    
    
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
    
    pop_agg3d[is.na(pop_agg3d)] = 0
    vc_agg3d[is.na(vc_agg3d)] = 0
    
    save(pop_agg3d, vc_agg3d, inc_v3d_agg, pop_moments_agg, file = "agg_pop_vc.RData", compress = "xz")
  } else {
    load("agg_pop_vc.RData")
  }
  
  print("agg_pop_vc")
  
  
  dim_survey = seroout$sero_studies
  dim_year = as.numeric(unique(pop1$year))
  dim_age = names(pop1 %>% select(-c(adm0_adm1, year, adm0)))

  #########################################################################################################
  ### LOAD PREVIOUS ###
  #########################################################################################################
  old_foi_est <- read.csv("multi_model_MCMC_chain_20190903/sero_chain_46_output5.csv", stringsAsFactors = FALSE)
  
  old_foi_est <- colMeans(old_foi_est[, 2:(ncol(old_foi_est)-4)])
  #########################################################################################################
  ### CREATE STATUS ###
  #########################################################################################################
  
  foi_const_surv = c(0,1e-6,0,0,0,0,rep(0,seroout$no_sero_surveys-6))
  
  out_p = create_pop_at_survey(pop_agg3d,
                               dim_survey = seroout$sero_studies,
                               dim_year = unique(pop1$year)) 
  p_at_survey = out_p$p_at_survey_3d
  P_tot_survey = out_p$P_tot_survey_2d
  
  
  # StartParamFoi=read.csv("sero_hpd_2019-08-21.csv",
  #                        header = TRUE)

  ## initialising parameters for MCMC
  parnames =  c("vac_eff",
                paste("Foi", seroout$sero_studies, sep = "_"),
                paste("R0", seroout$sero_studies, sep = "_"),
                paste("vc_factor_CMRs", sep = "_")
  )
  pars_ini = rep(NA,length(parnames))
  names(pars_ini) = parnames
  
  ## filling the pars_in vector
  # vaccine efficacy 
  pars_ini[1] = as.numeric(old_foi_est[grep("vac_eff", names(old_foi_est))])  # log(0.975) 
  
  # foi for each survey 
  pars_ini[grep("Foi", parnames)] = as.numeric(old_foi_est[grep("Foi", names(old_foi_est))]) #rep(log(0.01), seroout$no_sero_surveys)  
  
  # R0 for each survey KATY IS LOG TRANSFORMING THESE
  #this is the max post prob values of R0 from trial run
  pars_ini[grep("R0", parnames)] = rep(0.1, seroout$no_sero_surveys) 
  
  #vc.factor.CMRs 
  pars_ini[length(pars_ini)]= as.numeric(old_foi_est[grep("vc_factor", names(old_foi_est))]) #
  
  
  ## declare vector to identify different parameter types: vacc eff=1, Foi/R0=3, vc.factor.CMRs =4
  parameter_type = c(1,rep(3,2*seroout$no_sero_surveys), 4) # THIS NOW INDEXES THE DIFFERENT PARAMETER TYPES
  
  ## initial model 
  model_type = "Foi" #"R0" # 
  
  print("set pars_ini")
  #########################################################################################################
  ### LOAD PSEUDO PRIOR DISTRIBUTIONS ###
  #########################################################################################################
  # these don't matter if we are only fitting one model
  FOI_posterior_distributions <-
    read.csv(paste0("Z:/MultiModelInference/Foi/posterior_distributions_norm.csv"),
             stringsAsFactors = FALSE, 
             header = TRUE)
  FOI_posterior_distributions = cbind(FOI_posterior_distributions, model_type = "Foi")
  
  R0_posterior_distributions <- 
    read.csv(paste0("Z:/MultiModelInference/R0/posterior_distributions_norm_tweak.csv"), 
             stringsAsFactors = FALSE, 
             header = TRUE)
  R0_posterior_distributions = cbind(R0_posterior_distributions, model_type = "R0")
  
  posterior_distributions = rbind(FOI_posterior_distributions, R0_posterior_distributions)
  
  #tweak posterior distributions 
  scale_factor =  0.3 #setup_pseudoprior(pars_ini, posterior_distributions, "R0")  
  # 
  posterior_distributions$sd[posterior_distributions$model_type == "R0"] = 
  scale_factor * posterior_distributions$sd[posterior_distributions$model_type == "R0"]
  #########################################################################################################
  ### MCMC ###
  #########################################################################################################
  print("scaled posterior distributions")
  
  set.seed(run_id)
  
  ign = NA
  
  prob_Foi = log(1) 
  
  #create a directory to save the output in
  name_dir = paste0("new_prior_multi_model_MCMC_chain", "_", 
                    format(Sys.time(),"%Y%m%d"))
  dir.create(name_dir,  showWarnings = TRUE)
  
  Niter = 5e5
  
  print("ready to start estimation")
  
  Sero_Gibbs_MCMC(pars_ini = pars_ini,
                  seroout = seroout,
                  foi_const_surv = foi_const_surv,
                  vc_agg3d = vc_agg3d,
                  pop_agg3d = pop_agg3d,
                  pop_moments_agg=pop_moments_agg,
                  dim_year = dim_year,
                  dim_age = dim_age,
                  p_at_survey = p_at_survey,
                  P_tot_survey = P_tot_survey,
                  inc_v3d_agg = inc_v3d_agg,
                  model_type = model_type,
                  parameter_type = parameter_type,
                  prob_Foi = prob_Foi,
                  posterior_distributions = posterior_distributions,
                  ign = ign,
                  name_dir = name_dir,
                  Niter = Niter,
                  run_id = run_id)
  
}
