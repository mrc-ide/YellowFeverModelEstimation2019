
library(parallelsugar)

source("GLM_estimation.R")

model_var <- 1:20

mclapply(model_var, run_estimation, mc.cores = 4)
