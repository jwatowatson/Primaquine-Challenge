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
mod_master_pop = stan_model(file = 'Stan models/RBC_model_master_pop.stan')


## load trial data
load('RBC_model_data.RData')


# number of chains & iterations
nChains = 4
nIter = 4000
nthin = 8


####################### - Fit RBC model to data - #############################
stan_data = make_stan_dataset(my_data = PQdat)


if(i==1){
  mod_fit_pop1 = sampling(object = mod_master_pop, data = stan_data,
                          iter = nIter, chain = nChains, thin= nthin,seed=i,
                          init = make_init_list(4))
  save(mod_fit_pop1, PQ_dat_all, stan_data, ID_map, file = 'Rout/pop_fit1.RData')
}
if(i==2){
  stan_data2 = stan_data
  stan_data2$id = stan_data$id2
  mod_fit_pop2 = sampling(object = mod_master_pop, data = stan_data2,
                          iter = nIter, chain = nChains, thin= nthin,seed=i,
                          init = make_init_list(4))
  save(mod_fit_pop2, stan_data2, ID_map, file = 'Rout/pop_fit2.RData')
}



