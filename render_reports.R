
# script to render all the rmd for the different model variants
# should only be run after checking has converged 
rm(params)
#function to render each report quickly
render_report = function(model_variant, 
                         filename = "GLM_MCMC_diagnostics", 
                         foi_sample = FALSE) {

  rm()
  
  if(filename !="transmission_intensity_estimates"){
  rmarkdown::render(
    paste0(filename, ".Rmd"), 
    params = list(model_variant = model_variant),
    output_file = paste0(filename,  "-", model_variant, ".pdf")
  )
  } else {
    rmarkdown::render(
      paste0(filename, ".Rmd"), 
      params = list(model_variant = model_variant, foi_sample = foi_sample),
      output_file = paste0(filename,  "-", model_variant, ".pdf")
    )
  }
}


#-----------------------------------------------------------------------------
# glm
lapply(1:20, FUN = function(x) render_report(x, "GLM_MCMC_diagnostics"))

#-----------------------------------------------------------------------------
# transmission
lapply(1:20, FUN = function(x) render_report(x, "transmission_intensity_estimates",
                                            foi_sample = TRUE))

#-----------------------------------------------------------------------------
# burden
lapply(1:20, FUN = function(x) render_report(x, "burden_projections"))

#-----------------------------------------------------------------------------
# # compare with old
# lapply(1:20, FUN = function(x) render_report(x, "compare_with_old_foi"))
# 
# #-----------------------------------------------------------------------------
# # compare S America
# lapply(1:20, FUN = function(x) render_report(x, "Compare_SAmerica_burden"))


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
rmarkdown::render(
  paste0("GLM_comparison", ".Rmd"), 
  output_file = paste0("GLM_comparison", ".pdf")
)

rmarkdown::render(
  paste0("Foi_proportional_samples", ".Rmd"), 
  output_file = paste0("Foi_proportional_samples", ".pdf")
)

rmarkdown::render(
  paste0("FOI_comparison", ".Rmd"), 
  output_file = paste0("FOI_comparison", ".pdf")
)

for(y in c(1995, 2005, 2013, 2018)){
  rmarkdown::render(
    paste0("burden_proportional_samples", ".Rmd"), 
    params = list(year = y),
    output_file = paste0("burden_proportional_samples_",y , ".pdf")
  )
}

rmarkdown::render(
  paste0("burden_compare", ".Rmd"), 
  output_file = paste0("burden_compare", ".pdf")
)

#-----------------------------------------------------------------------------
for(y in c(2012, 2013, 2018)){
  rmarkdown::render(
    paste0("vaccine_impact_projections", ".Rmd"), 
    params = list(year = y),
    output_file = paste0("vaccine_impact_projections_",y , ".pdf")
  )
}
