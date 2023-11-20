fun_stoch = function(run_id){
  
  param_samples = FOIestimates[run_id,]
  
  for(c in 1:length(coverage_types)){
    
    if(!file.exists(paste0(outdir,"/","stochastic-burden-estimates.",
                           touchstone_version, "_YF_IC-Garske_", 
                           coverage_types[c], "_", run_id, ".csv"))){
      
      GAVI_vac <- GAVI_vac_list[[c]]
      
      if(anyNA(as.numeric(GAVI_vac$coverage))){
        GAVI_vac[is.na(as.numeric(GAVI_vac$coverage)),"coverage"] <- 0
      }
      
      GAVI_vac %<>% mutate(coverage = as.numeric(coverage))
      
      
      ### run stochastic ###----------------------------------------------------------------------------
      output_df = run_burden_for_template(historic_dat=NULL, #historic now incorporated into Montagu
                                          GAVI_vac,
                                          GAVI_switch = "preventive", 
                                          param_samples,
                                          vac_eff = vac_eff[run_id],
                                          pop_all,
                                          P_severe[run_id],
                                          P_severeDeath[run_id],
                                          life_exp_GAVI,
                                          template,
                                          round = TRUE,
                                          run_id = run_id)
      
      
      # checks
      names(output_df)
      names(template)
      nrow(template)==nrow(output_df)
      
      
      write.csv(output_df, 
                paste0(outdir,"/","stochastic-burden-estimates.",
                       touchstone_version, "_YF_IC-Garske_", 
                       coverage_types[c], "_", 
                       run_id, ".csv"), 
                row.names = FALSE)
      
    }
  }
}