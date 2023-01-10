args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file



## Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## for parallelisation
cores = detectCores()

source('master_functions.R')

## load Stan model
mod = stan_model(file = 'Stan models/RBC_model_master_pop_CVpred.stan')


load('RBC_model_data.RData')

# set number of subjects
IDs = unique(PQ_dat_all$ID)
# number of chains & iterations

nChains = 4
nIter = 4000
nthin = 8

# remove patient i and add to predict
if(i < max(PQ_dat_all$ID)){
  PQ_dat_pred = PQ_dat_all[PQ_dat_all$ID == i, ]
  PQ_dat_fit = PQ_dat_all[PQ_dat_all$ID != i, ]
  
  Tmax_pred = max(PQ_dat_pred$Study_Day) - min(PQ_dat_pred$Study_Day) + 1 
  drug_regimen_pred = array(0, dim = Tmax_pred) # initialise with all zeros
  
  timepoints = PQ_dat_pred$Study_Day - min(PQ_dat_pred$Study_Day) + 1
  drug_regimen_pred[timepoints] = PQ_dat_pred$dosemgkg
  
}
####################### - Fit RBC model to data - #############################
stan_data = data_model(my_data = PQ_dat_fit,
                       N_pred = length(drug_regimen_pred),
                       drug_regimen_pred = drug_regimen_pred)

out = sampling(object = mod_master_pop, data = stan_data,
               iter = nIter, chain = nChains, thin= nthin,seed=i,
               init = make_init_list(4))
save(out, file = paste0('Rout/CV_fit',i,'.RData'))




