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
mod_master_pop = stan_model(file = 'Stan_models/RBC_single_individ_mechanistic_G6PD.stan')

# number of chains & iterations
nChains = 4
nIter = 4000
nthin = nChains

load('Rout/stan_data_list.RData')
if(job_i > length(dat_stan_list)) stop('no job to do')

dat_stan_list[[job_i]]$logit_MAX_EFFECT_prior_mean = -4.5
dat_stan_list[[job_i]]$logit_MAX_EFFECT_prior_sigma = 0.5
out = sampling(object = mod_master_pop, 
               data = dat_stan_list[[job_i]],
               iter = nIter, chain = nChains, thin = nthin,
               seed = job_i,
               init = make_init_list(nChains),
               pars = c('L_Omega','L_Omega_ic'), include = FALSE)

traceplot(out, pars = c('diff_alpha', 'delta_alpha', 
                        'logit_G6PD_delta_day','logit_MAX_EFFECT','logit_G6PD_threshold',
                        'Hb_star', 'h','beta','log_k' ))

thetas=extract(out, pars='Y_pred')$Y_pred

par(mfrow=c(1,2), las=1)
plot(colMeans(thetas[,1,]), type='l',lwd=3)
points(dat_stan_list[[job_i]]$t_sim_hemocue_Hb,
       dat_stan_list[[job_i]]$Hb_Haemocue, col='red',pch=15)

plot(colMeans(thetas[,2,]), type='l',lwd=3)
points(dat_stan_list[[job_i]]$t_sim_retic,
       dat_stan_list[[job_i]]$Retic_data, col='red',pch=15)

save(out, file = paste0('Rout/job_',job_i,'.RData'))

writeLines(sprintf('Finished job %s', job_i))
