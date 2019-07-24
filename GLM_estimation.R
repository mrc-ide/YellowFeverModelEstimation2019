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

dat %<>% filter(!is.na(precip_mean))

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# FIT MODEL #

object_glm = YFestimation::fit_glm(dat = dat, 
                                   depi = match("cases_or_outbreaks", names(dat)), 
                                   models = "cases_or_outbreaks~surv.qual.adm0+adm05+logpop+temp_max+precip_mean+altitude" )  

beta0 = object_glm[[1]]
x = object_glm[[2]]
y = object_glm[[3]]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# INITIAL PARAMETERS #

#pars_ini = beta0


pars_ini = as.numeric(read.csv("GLM_parameters_2019-07-23.csv", stringsAsFactors = TRUE))
pars_ini = c(pars_ini, beta0[length(beta0)])
names(pars_ini) = paste0("log.",names(beta0))

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
name_dir = paste0("GLM_MCMC_chain", "_", format(Sys.time(),"%Y%m%d"))
dir.create(name_dir)

Niter = 1e6

# MCMC #
YFestimation::GLM_MCMC(Niter, name_dir, pars_ini, x, y, plot_chain)
