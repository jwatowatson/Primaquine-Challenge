args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
job_i = as.numeric(args[6])
print(paste0("job(i) = ", job_i)) # this will print out in the *.o file


## Load required packages
packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
lapply(packs, require, character.only = TRUE)

options(mc.cores = parallel::detectCores())

source('master_functions.R')

## load Stan model
mod_master_pop = stan_model(file = 'Stan_models/RBC_model_mechanistic_G6PD.stan')

# number of chains & iterations
nChains = 1
nIter = 1000
nthin = nChains

load('Rout/stan_data_list.RData')
if(job_i > length(dat_stan_list)) stop('no job to do')

dat_stan_list[[job_i]]$log_MAX_EFFECT_prior_mean = -2.5
dat_stan_list[[job_i]]$log_MAX_EFFECT_prior_sigma = 0.5
dat_stan_list[[job_i]]$log_beta_mean = -1.5
dat_stan_list[[job_i]]$log_k_sigma = 0.25
dat_stan_list[[1]]$T_RBC_max=140

out = sampling(object = mod_master_pop, 
               data = dat_stan_list[[job_i]],
               iter = nIter, chain = nChains, thin = nthin,
               seed = job_i,verbose=T, 
               # init = make_init_list(nChains),
               control=list(max_treedepth=10),
               pars = c('L_Omega'), include = FALSE)

shinystan::launch_shinystan(out)
sampler_params <- get_sampler_params(out, inc_warmup = F)
summary(do.call(rbind, sampler_params), digits = 2)
pairs(out, pars = c("log_MAX_EFFECT", "log_G6PD_decay_rate"), las = 1)

pairs(out, pars = c('diff_alpha', 'delta_alpha', 
                    'log_G6PD_decay_rate','log_MAX_EFFECT',
                    'sigma_death',
                    'Hb_star', 'h','log_beta','log_k'), las = 1)
traceplot(out, pars = c('diff_alpha', 'delta_alpha', 
                        'log_G6PD_decay_rate','log_MAX_EFFECT',
                        'sigma_death',
                        'Hb_star', 'h','log_beta','log_k',
                        'sigma_CBC',
                        'sigma_haemocue',
                        # 'sigma_retic',
                        'CBC_correction'),inc_warmup=F)

summary(out, pars = c('diff_alpha', 'delta_alpha', 
                        'log_G6PD_decay_rate','log_MAX_EFFECT',
                        'sigma_death',#'mu_death',
                        'Hb_star', 'h','log_beta','log_k',
                        'sigma_CBC',
                        'sigma_haemocue',
                        'sigma_retic',
                        'CBC_correction'))$summary
thetas=extract(out, pars='Y_hat')$Y_hat
plot(dat_stan_list[[job_i]]$t_sim_hemocue_Hb,
     colMeans(thetas[,1,]))
plot(dat_stan_list[[job_i]]$t_sim_hemocue_Hb,
     dat_stan_list[[job_i]]$Hb_Haemocue-colMeans(thetas[,1,]))

ind_retic = get_ind_retic(dat_stan_list[[job_i]])
plot(dat_stan_list[[job_i]]$t_sim_retic,
     colMeans(thetas[,2,ind_retic]))

plot( colMeans(thetas[,3,ind_retic]))
plot( colMeans(thetas[,4,ind_retic]))

thetas=extract(out, pars='Y_pred')$Y_pred

par(mfrow=c(1,2), las=1)
plot(colMeans(thetas[,1,]), type='l',lwd=3, ylim =c(8,16))
points(dat_stan_list[[job_i]]$t_sim_hemocue_Hb,
       dat_stan_list[[job_i]]$Hb_Haemocue, col='red',pch=15)

plot(colMeans(thetas[,2,]), type='l',lwd=3, ylim=c(0,15))
points(dat_stan_list[[job_i]]$t_sim_retic,
       dat_stan_list[[job_i]]$Retic_data, col='red',pch=15)
plot(colMeans(thetas[,3,]), type='l',lwd=3)
plot(colMeans(thetas[,4,]), type='l',lwd=3)

save(out, file = paste0('Rout/job_',job_i,'.RData'))

writeLines(sprintf('Finished job %s', job_i))
