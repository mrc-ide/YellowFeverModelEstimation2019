
library(parallelsugar)

source("GLM_estimation.R")

run_id = 1:4

mclapply(run_id, run_estimation, mc.cores = 4)
