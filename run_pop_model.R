args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file



## Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

#### set seed (mm/YY) ####
set.seed(0522)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## for parallelisation
cores = detectCores()

source('master_functions.R')

## load Stan model
mod_master_pop = stan_model(file = 'Stan models/RBC_model_master_pop_simplify.stan')

# number of chains & iterations
nChains = 1
nIter = 1000
nthin = nChains

####################### - Fit RBC model to data - #############################
## load trial data
load('RBC_model_data.RData')
stan_data = make_stan_dataset(my_data = PQdat)

mod_fit_pop1 = sampling(object = mod_master_pop, 
                        data = stan_data,
                        iter = nIter, chain = nChains, thin = nthin,
                        seed = 1234,
                        init = make_init_list(nChains),
                        pars = 'L_Omega', include = FALSE)
save(mod_fit_pop1, PQdat, stan_data, file = 'Rout/pop_fit1.RData')



traceplot(mod_fit_pop1, pars=pars)

preds = rstan::extract(mod_fit_pop1, pars='Y_pred')$Y_pred


