args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
job_i = as.numeric(args[6])
print(paste0("job(i) = ", job_i)) # this will print out in the *.o file

job_i=5

## Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

options(mc.cores = parallel::detectCores())

source('master_functions.R')

## load Stan model
mod_master_pop = stan_model(file = 'Stan_models/RBC_model_mechanistic_G6PD.stan')

# number of chains & iterations
nChains = 4
nIter = 2000
nthin = nChains

load('Rout/stan_data_list.RData')
if(job_i > length(dat_stan_list)) stop('no job to do')

dat_stan_list[[job_i]]$logit_MAX_EFFECT_prior_mean = -2
dat_stan_list[[job_i]]$logit_MAX_EFFECT_prior_sigma = 0.1
out = sampling(object = mod_master_pop, 
               data = dat_stan_list[[job_i]],
               iter = nIter, chain = nChains, thin = nthin,
               seed = job_i,
               init = make_init_list(nChains),
               pars = c('L_Omega','L_Omega_ic'), include = FALSE)

traceplot(out, pars = c('diff_alpha', 'delta_alpha', 'logit_G6PD_delta_day',
                        'sigma_CBC','sigma_haemocue','sigma_retic',
                        'CBC_correction','Hb_star','logit_MAX_EFFECT',
                        'h','beta','log_k', 'logit_G6PD_threshold'))

thetas=extract(out, pars='Y_pred')$Y_pred
plot(colMeans(thetas[,1,]))
plot(colMeans(thetas[,2,]))

save(out, file = paste0('Rout/job_',job_i,'.RData'))

writeLines(sprintf('Finished job %s', job_i))
