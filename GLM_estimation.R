#mcmc for estimating glm

run_estimation = function(run_id =1){
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # LIBRARIES FOR PACKAGES USED #
  
  library(dplyr)
  library(readr)
  library(reshape2)
  library(mvtnorm)
  library(YFestimation)
  library(KsetupR)
  library(R.utils)
  library(magrittr)
  library(truncdist)
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # NEW FUNCTIONS #
  sourceDirectory("Functions", modifiedOnly = FALSE)
  
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # LOAD ENVIRONMENTAL DATA #
  
  Env_Table_path = "../Data/Environment/global_dat"
  
  filename = get_latest_file(path = Env_Table_path, pattern = "dat_wes")
  
  dat = read.csv(filename, stringsAsFactors = FALSE)
  
  # remove families of NHP that are not to be included
  model_form_whole = read.csv("Model_form_a.csv", stringsAsFactors = FALSE)$x
  covar = unlist(strsplit(unlist(strsplit(model_form_whole, "\\+")), "\\~"))
  
  dat = dat[, -which(!names(dat) %in% covar & grepl("family", names(dat)))] # remove family which are not covariates

  #try fitting just to high/moderate
  #dat %<>% filter(risk %in% c("high", "moderate"))
  
  dat %<>% dplyr::select(-year)
  dat %<>% filter( is.finite(MIR.max))
  dat %<>% filter(!is.na(temp_mean))
  #dat %<>% tidyr::drop_na()
  dat %<>% mutate(continent = ifelse(adm0 %in% c("CIV", "STP"),
                                     "Africa", continent))
  
  # make extra elements and normalise
  dat = adjust_env_dat(dat)
  
  
  
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # FIT MODEL #
  covar = covar[grep("cases_or_outbreaks|family", covar, invert = TRUE)]
  
  model_form = paste0("cases_or_outbreaks~", paste(covar, collapse = "+"),"+aggregate_family" ,  "+adm05", "+surv.qual.adm0")# 

  
  
  object_glm = YFestimation::fit_glm(dat = dat, 
                                     depi = match("cases_or_outbreaks", names(dat)), 
                                     models = model_form) #"+predicted_surv_qual"))  
  
  beta0 = object_glm[[1]]
  x = object_glm[[2]]
  y = object_glm[[3]]

  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # INITIAL PARAMETERS #
  
  pars_ini = beta0
  pars_ini[grep("family", names(pars_ini))] = abs(pars_ini[grep("family", names(pars_ini))])

  pars_ini[grep("low_risk", names(pars_ini))] = 1
  names(pars_ini) = paste0("log.", names(beta0))
  
  # parm_in = read.csv("glm_param.csv", stringsAsFactors = TRUE)
  # parm_in = parm_in[1:(length(parm_in)-2)]
  # pars_ini = as.numeric(parm_in)
  # names(pars_ini) = paste0("log.", names(parm_in))

  pars_ini = pars_ini*0
  
  pars_ini[is.na(pars_ini)] = 0
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # SET SEED #
  index_code = as.numeric(commandArgs(TRUE))
  if(identical(index_code, numeric(0))) index_code=1
  
  t= as.numeric(Sys.time())
  seed= (t - floor(t)) * 1e8 
  set.seed(seed)
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # MCMC #
  
  # set-up #
  plot_chain = FALSE
  
  # create a directory to save the output in 

  name_dir = paste0("GLM_MCMC_chain", "_", format(Sys.time(),"%Y%m%d"), "_a")

  dir.create(name_dir)
  
  Niter = 1e6
  
  # MCMC #
  GLM_MCMC(Niter, name_dir, pars_ini, x, y, plot_chain, run_id)
  
}
