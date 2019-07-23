#mcmc for estimating glm

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

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# NEW FUNCTIONS #
sourceDirectory("Functions", modifiedOnly = FALSE)


# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD ENVIRONMENTAL DATA #

Env_Table_path = "../Data/Environment/global_dat"

filename = get_latest_file(path = Env_Table_path, pattern = "dat")

dat = read.csv(filename, stringsAsFactors = FALSE)

dat = adjust_env_dat(dat)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# FIT MODEL #

object_glm = YFestimation::fit_glm(dat = dat, 
                                   depi = match("cases_or_outbreaks", names(dat)), 
                                   models = "cases_or_outbreaks~surv.qual.adm0+adm05+logpop+temp_max+precip_mean" )  

beta0 = object_glm[[1]]
x = object_glm[[2]]
y = object_glm[[3]]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# INITIAL PARAMETERS #