
library(mcmcplots)
library(ggmcmc)
library(gridExtra)
library(maptools)
library(readr)
library(fields)
library(Hmisc)
library(dplyr)
library(magrittr)
library(R.utils)

library(YFestimation)
library(snapalette)
library(KsetupR)


sourceDirectory("Functions", modifiedOnly = FALSE)

snapal <- "Pop"

filepath_to_use <- "multi_model_MCMC_chain_20191129"

mcmc_out <- get_chains(filepath_to_use,
                       burnin = 1, thin = 100)

filepath_new_prior <- "new_prior_multi_model_MCMC_chain_20210128"

mcmc_out_new_prior <- get_chains(filepath_new_prior,
                       burnin = 1, thin = 100)

mcmc_out_new_prior %<>% mutate(prior_type = "rate = 0.1")
mcmc_out %<>% mutate(prior_type = "rate = 0.001")

mcmc_all <- bind_rows(mcmc_out, mcmc_out_new_prior)

mcmc_all  %<>% select(-starts_with("R0"))

mcmc_all_long <- mcmc_all %>%
  mutate(sample = 1:nrow(mcmc_all)) %>%
  pivot_longer(names_to = "parameter", values_to = "estimate", -c(sample, prior_type))


mcmc_all_long %>%
  filter(grepl("Foi", parameter)) %>%
  mutate(parameter = gsub("Foi_", "", parameter)) %>%
  ggplot()+
  geom_boxplot()+
  aes(x = parameter, 
      y = exp(estimate),
      fill = prior_type)+
  scale_y_log10()+
  theme_bw()+
  labs(x = "Serology study site", y = "Force of infection estimate", fill = "FOI prior")+
  scale_fill_manual(values = snapalette(snapal,5)[c(2,5)])+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

