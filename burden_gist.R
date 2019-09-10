# burden gist

library(YFburden)
library(dplyr)
library(magrittr)
library(tidyr)
library(snapalette)
library(KsetupR)
library(R.utils)
library(YFestimation)

sourceDirectory("Functions", modifiedOnly = FALSE)

# ---------------------------------------------------------------------------
### POPULATION DATA ###

all_res_pop_3d = get_pop_data_3d() 

pop1 = all_res_pop_3d$pop1                                      
P_tot = all_res_pop_3d$P_tot                                   #total populations for each adm and year
p_prop = all_res_pop_3d$p_prop                                   #proportions of population

pop1 %<>% rename(adm0_adm1 = adm1)

pop1 %<>% arrange(year) %>% select(-c(country_code, country))

names(pop1)[3:ncol(pop1)] = paste0("a", 0:100)

pop1 %<>% filter(year<2051)
p_prop %<>% filter(year<2051)
P_tot %<>% filter(year<2051)

# ---------------------------------------------------------------------------
# VACCINATION DATA #
vc2d = readRDS("vacc_1940_1950.RDS")
vc2d = as.data.frame(vc2d)

# expand vc2d to match pop1
vc2d %<>% mutate(year = as.numeric(as.character(year)))
vc_tmp = pop1 %>% filter(!adm0_adm1 %in% vc2d$adm0_adm1)
vc_tmp[, grep("adm0_adm1|year", names(vc_tmp), invert = TRUE)] = 0

vc2d %<>% bind_rows(vc_tmp)
vc2d %<>% unique()
vc2d %<>% filter(!is.na(adm0_adm1))


# ---------------------------------------------------------------------------
# FOI #
foi = read.csv( "FOI_samples.csv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
fun.calcBurden = function(fois,year,adm1 = NA) { # Tini's burden function
  # using vc.full and pop.full as global variables.
  if(is.na(adm1[1])) adm1 = pop.full$adm1[pop.full$year==year]
  mm = match(adm1,pop.full$adm1[pop.full$year==year])
  
  foi.mat = matrix(fois,nrow=length(fois),ncol=101,byrow=F)
  a.mat = matrix(0:100,nrow=length(fois),ncol=101,byrow=T)
  burden = fois*rowSums(exp(-foi.mat*a.mat)*(1-vc.full[vc.full$year==year,-(1:3)])*
                          pop.full[pop.full$year==year,-(1:3)], na.rm=T)
  return(burden)
}
