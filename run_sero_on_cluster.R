#### run on cluster ####

# drat:::addRepo("mrc-ide")
# install.packages("didehpc")
library(didehpc)

didehpc_config()

context::context_log_start()

root = "contexts"

packages <- list( attached = c("dplyr", "readr", "reshape2", "mvtnorm", "R.utils",  "YFestimation", "magrittr"))

sources <- c("Functions/custom_estimation.R", 
             "Functions/data_launch.R", 
             "Functions/diagnostics.R", 
             "sero_estimation.R")

ctx <- context::context_save(root, packages = packages, sources = sources,
                             package_sources = provisionr::package_sources(local = c( "Z:/YFestimation", "Z:/Burden_package/YFburden") ))


obj <- didehpc::queue_didehpc(ctx)


# run for a number

run_id = 10:14
grp <- obj$lapply(run_id, run_estimation)