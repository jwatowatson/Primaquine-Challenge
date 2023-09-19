#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

# Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("master_functions.R")

# Load the Stan model
mod_master_pop <- stan_model(
  file = "Stan models/RBC_model_master_pop_free_weights.stan"
)

# Number of chains & iterations
chains <- 4
iters <- 1000
nthin <- chains

# Fit the model to data
load("Data/RBC_model_data.RData")
stan_data <- make_stan_dataset(my_data = PQdat)

mod_fit_pop_free_wt <- sampling(
  object = mod_master_pop,
  data = stan_data,
  iter = iters,
  chain = chains,
  thin = nthin,
  seed = 1234,
  init = make_init_list(chains),
  pars = c("L_Omega", "L_Omega_ic"),
  include = FALSE
)

# Save the results
results_file <- "Rout/pop_fit_free_weights.RData"
save(mod_fit_pop_free_wt, PQdat, stan_data, file = results_file)
