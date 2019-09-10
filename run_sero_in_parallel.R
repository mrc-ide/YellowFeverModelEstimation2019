
library(parallelsugar)

source("sero_estimation.R")

run_id = 2:5

mclapply(run_id, run_estimation, mc.cores = 4)
