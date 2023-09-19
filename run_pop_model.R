args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
job_i = as.numeric(args[6])
print(paste0("job(i) = ", job_i)) # this will print out in the *.o file



## Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

options(mc.cores = parallel::detectCores())

source('master_functions.R')

## load Stan model
mod_master_pop = stan_model(file = 'Stan_models/RBC_model_master_pop_free_weights.stan')

# number of chains & iterations
nChains = 4
nIter = 1000
nthin = nChains

load('Rout/stan_data_list.RData')
if(job_i > length(dat_stan_list)) stop('no job to do')

dat_stan_list[[job_i]]$K_weights = 14
dat_stan_list[[job_i]]$prior_weights = c(rep(10,7), rep(5, 3), rep(1,4))

out = sampling(object = mod_master_pop, 
               data = dat_stan_list[[job_i]],
               iter = nIter, chain = nChains, thin = nthin,
               seed = job_i,
               init = make_init_list(nChains),
               pars = c('L_Omega','L_Omega_ic'), include = FALSE)
save(out, file = paste0('Rout/job_',job_i,'.RData'))

writeLines(sprintf('Finished job %s', job_i))
