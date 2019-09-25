###########################################################
### FUNCTIONs for temp suitability ###
###########################################################
briere = function(Temp, parm){
  out = parm[3]*Temp*(Temp-parm[1])*(parm[2]-Temp)^0.5
  return(out)
}

quad = function(Temp, parm){
  out = -parm[3]*(Temp-parm[1])*(parm[2] -Temp)
  return(out)
}
###########################################################

###########################################################
### FUNCTION to calculate temperature suitability index ###
###########################################################

temp_suitability = function(Temp, param){
  
  #bite rate
  a = briere(Temp, 
             as.numeric(param[ grep("^a_", names(param))] )) # new estimate medians
  
  # mu = m mortality
  lf = quad(Temp, 
            as.numeric(param[ grep("^mu_", names(param))]))
  mu = 1/ ifelse( lf<=0, 1, lf ) #guard against negatives and zeros
  
  
  # PDR = parasite development rate = 1/EIP  
  PDR = briere(Temp, 
               as.numeric(param[ grep("^PDR_", names(param))])) #
  PDR[PDR<0]=0
  
  a[is.na(a)] = PDR[is.na(PDR)] = 0
  
  # suitability
  Z = (a^2 * exp(-mu / PDR) ) / mu 
  
  return(Z)
}