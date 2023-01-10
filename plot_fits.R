library(gtools)
library(rstan)
source('master_functions.R')
load('Rout/pop_fit1.RData')
load('Rout/pop_fit2.RData')

pars = c("Hb_steady_state", "rbc_lifespan_steady_state",
         "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
         "logit_max_drug_effect",
         "ed50", "emax_k",'mean_delay','logit_max_drug_effect',
         'theta_rand')

traceplot(mod_fit_pop1, 
          pars=c("Hb_steady_state", "rbc_lifespan_steady_state",
                 "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
                 "logit_max_drug_effect",'rbc_lifespan_steady_state',
                 "ed50", "emax_k",'mean_delay','log_k',
                 'sigma_delay','logit_max_drug_effect',
                 'sigma_CBC', 'sigma_manual'))

traceplot(mod_fit_pop2, 
          pars=c("Hb_steady_state", "rbc_lifespan_steady_state",
                 "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
                 "logit_max_drug_effect",
                 "ed50", "emax_k",'mean_delay',
                 'sigma_delay','logit_max_drug_effect',
                 'sigma_CBC', 'sigma_manual'))


Y_hat1 = rstan::extract(mod_fit_pop1, pars='Y_hat')$Y_hat
Y_hat2 = rstan::extract(mod_fit_pop2, pars='Y_hat')$Y_hat

my_thetas1 = rstan::extract(mod_fit_pop1, pars=pars)
my_thetas2 = rstan::extract(mod_fit_pop2, pars=pars)

# pdf('model_fits_all.pdf',width = 10, height = 10)
par(mfrow=c(2,2),las=1,family='serif',cex.lab=1.3, cex.axis=1.3)

for(id in 1:stan_data$N){
  ind_id = stan_data$ind_start_regimen[id]:stan_data$ind_end_regimen[id]
  Y_hat_hb = Y_hat1[,1,ind_id]
  Y_hat_retic = Y_hat1[,2,ind_id]
  
  Y_hat_hb2 = Y_hat2[,1,ind_id]
  Y_hat_retic2 = Y_hat2[,2,ind_id]
  
  
  ind_haemocue = stan_data$t_sim_hemocue_Hb[stan_data$ind_start_Hb_hemocue[id]:stan_data$ind_end_Hb_hemocue[id]]
  ind_hb_CBC = stan_data$t_sim_CBC_Hb[stan_data$ind_start_Hb_CBC[id]:stan_data$ind_end_Hb_CBC[id]]
  ind_retic = stan_data$t_sim_retic[stan_data$ind_start_retic[id]:stan_data$ind_end_retic[id]]
  
  plot(1:length(ind_id), colMeans(Y_hat_hb), type='l', lwd=3,
       panel.first=grid(), ylim =c(8.3,16.3), xlab='Time (days)',
       ylab='Haemoglobin (g/dL)',main=paste0('Subject ',id))
  lines(1:length(ind_id), colMeans(Y_hat_hb2))
  dosing_ind = which(stan_data$drug_regimen[ind_id]>0)
  dosing_ind = c(head(dosing_ind,1), tail(dosing_ind,1)+1)
  polygon(c(dosing_ind, rev(dosing_ind)), c(0,0,20,20),col=adjustcolor('pink',.3),border = NA)
  
  polygon(c(1:length(ind_id), rev(1:length(ind_id))),
          c(apply(Y_hat_hb,2, quantile,probs=0.1), rev(apply(Y_hat_hb,2,quantile,probs=0.9))),
          col=adjustcolor('grey',.3),border = NA)
  lines(1:length(ind_id), colMeans(Y_hat_hb), lwd=3)
  points((1:length(ind_id))[ind_haemocue],
         stan_data$Hb_Haemocue[stan_data$ind_start_Hb_hemocue[id]:stan_data$ind_end_Hb_hemocue[id]], 
         pch =15,cex=1.3)
  points((1:length(ind_id))[ind_hb_CBC],
         stan_data$Hb_CBC[stan_data$ind_start_Hb_CBC[id]:stan_data$ind_end_Hb_CBC[id]], 
         pch =16, col='blue',cex=1.3)
  
  plot(1:length(ind_id), colMeans(Y_hat_retic), type='l', lwd=3,
       panel.first=grid(), ylim =c(0,18), xlab='Time (days)',
       ylab='Reticulocyte count (%)',main=paste0('Subject ',id))
  polygon(c(dosing_ind, rev(dosing_ind)), c(0,0,20,20),col=adjustcolor('pink',.3),border = NA)
  lines(1:length(ind_id), colMeans(Y_hat_retic), lwd=3)
  points((1:length(ind_id))[ind_retic],
         stan_data$Retic_data[stan_data$ind_start_retic[id]:stan_data$ind_end_retic[id]], 
         pch =15,cex=1.3)
}

par(mfrow=c(1,1))

xs= seq(0.01,1,length.out=200)
effect = array(dim=c(length(my_thetas1$ed50), length(xs)))
for(i in 1:length(xs)){
  effect[,i] = emax_f(x = xs[i], 
                      logit_max_drug_effect = my_thetas1$logit_max_drug_effect,
                      emax_k = my_thetas1$emax_k, 
                      ed50 = my_thetas1$ed50)
}
effect_mean = colMeans(effect)
plot(xs, 100*effect_mean,type='l',lwd=3, xlab='Primaquine mg/kg',
     ylab='Reduction in red cell lifespan (%)',panel.first=grid(),
     ylim = c(0,100))
polygon(c(xs, rev(xs)),
        100*c(apply(effect,2, quantile,probs=0.1), rev(apply(effect,2,quantile,probs=0.9))),
        col=adjustcolor('grey',.3),border = NA)

effect2 = array(dim=c(length(my_thetas2$ed50), length(xs)))
for(i in 1:length(xs)){
  effect2[,i] = emax_f(x = xs[i], 
                      logit_max_drug_effect = my_thetas2$logit_max_drug_effect,
                      emax_k = my_thetas2$emax_k, 
                      ed50 = my_thetas2$ed50)
}
effect_mean2 = colMeans(effect2)
lines(xs, effect_mean2*100)

effect_id_list = effect_id_list2 = list()
for(id in 1:stan_data$N){
  logit_max_drug_effect = my_thetas1$logit_max_drug_effect+my_thetas1$theta_rand[,id,6]
  emax_k = my_thetas1$emax_k*exp(my_thetas1$theta_rand[,id,7])
  ed50 = my_thetas1$ed50*exp(my_thetas1$theta_rand[,id,8])
  
  j = stan_data2$id[id]
  logit_max_drug_effect2 = my_thetas2$logit_max_drug_effect+my_thetas2$theta_rand[,j,6]
  emax_k2 = my_thetas2$emax_k*exp(my_thetas1$theta_rand[,j,7])
  ed502 = my_thetas2$ed50*exp(my_thetas1$theta_rand[,j,8])
  
  effect_id = effect_id2 = array(dim=c(length(my_thetas1$ed50), length(xs)))
  for(i in 1:length(xs)){
    effect_id[,i]= emax_f(x = xs[i], 
                          logit_max_drug_effect = logit_max_drug_effect,
                          emax_k = emax_k, ed50 = ed50)
    effect_id2[,i]= emax_f(x = xs[i], 
                          logit_max_drug_effect = logit_max_drug_effect2,
                          emax_k = emax_k2, ed50 = ed502)
  }
  effect_id_list[[id]]=colMeans(effect_id)
  effect_id_list2[[j]]=colMeans(effect_id2)
  
  lines(xs, 100*colMeans(effect_id))
}
lines(xs, 100*colMeans(effect),lwd=3,col='red')

ids_both = ID_map$ID_combined[duplicated(ID_map$ID_combined)]
par(mfrow=c(3,3))
for(id in ids_both){
  plot(xs, 100*colMeans(effect),type='l',lwd=3, xlab='Primaquine mg/kg',
       ylab='Reduction in red cell lifespan (%)',panel.first=grid(),
       ylim = c(0,60))
  is = which(ID_map$ID_combined==id)
   lines(xs, 100*effect_id_list[[is[1]]], lwd=2, lty=2,col='red')
   lines(xs, 100*effect_id_list[[is[2]]], lwd=2, lty=2,col='green')
  
  lines(xs, 100*effect_id_list2[[is[1]]], lwd=2, lty=2,col='blue')
}

theta_rand_mean = apply(my_thetas$theta_rand, 2:3, mean)
ID_map$Hb_star= mean(my_thetas$Hb_steady_state)+theta_rand_mean[,1]
plot(ID_map$Hb_star, ID_map$Hb_star, type='n',panel.first=grid())
for(id in ids_both) {
  is = which(ID_map$ID_combined==id)
  points(ID_map$Hb_star[is[1]], ID_map$Hb_star[is[2]])
}

dev.off()


xx = my_thetas$theta_rand[,,9]
for(i in 1:ncol(xx)){
  xx[,i] = xx[,i]+my_thetas$rbc_lifespan_steady_state
}
hist(colMeans(xx),breaks = seq(30,90,by=5))
