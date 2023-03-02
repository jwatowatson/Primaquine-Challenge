library(gtools)
library(rstan)
source('master_functions.R')
load('Rout/pop_fit1.RData')

plot_model_fit = function(mod_fit, stan_data) {
  pars = c('mean_delay','sigma_delay','sigma_CBC','sigma_haemocue','logit_alpha',
           'h','beta','log_k','Hb_star','T_E_star','sigmasq_u',
           'alpha_diff1','alpha_diff2','alpha_delta1','alpha_delta2' )
  traceplot(mod_fit, pars=pars)
  preds = rstan::extract(mod_fit, pars = 'Y_pred')$Y_pred
  
  summary_pred = apply(preds, 2:3, quantile, probs = c(0.1, .5, .9))
  par(mfrow=c(2,2),las=1, family='serif')
  plot(summary_pred[2,1, ],type='l', ylim = range(preds[,1,]), lwd=3,
       panel.first=grid(), xlab='Days', ylab='Haemoglobin (g/dL)')
  lines(summary_pred[3,1, ])
  lines(summary_pred[1,1, ])
  
  plot(summary_pred[2,2, ],type='l', ylim = c(0,20), lwd=3,
       panel.first=grid(), xlab='Days', ylab='Reticulocytes (%)')
  lines(summary_pred[3,2, ])
  lines(summary_pred[1,2, ])
  
  pred_falls = apply(preds[,1,], 1, function(x) x[1] - min(x))
  hist(pred_falls, xlab='Haemoglobin fall (g/dL)', ylab='', 
       breaks=20, yaxt='n',main='')
  
  plot(preds[which.max(pred_falls),1,], type='l', lwd=3,
       panel.first=grid(), xlab='Days',ylab='Haemoglobin (g/dL)')
  
  
  # Residuals plot
  Y_hat = rstan::extract(mod_fit, pars = 'Y_hat')$Y_hat 
  ## CBC
  Hb_CBC_observed = stan_data$Hb_CBC
  ind_cbc = get_ind_cbc(stan_data)
  Hb_CBC_mean_predicted = colMeans(Y_hat[,1,ind_cbc]) + 
    mean(extract(mod_fit, pars='CBC_correction')$CBC_correction)
  hist(Hb_CBC_observed-Hb_CBC_mean_predicted,
       xlab='Residuals (CBC haemoglobin)',main='')
  boxplot( (Hb_CBC_observed-Hb_CBC_mean_predicted) ~ stan_data$t_sim_CBC_Hb)
  abline(h=0, lwd=2, col='red')
  
  ## Haemocue
  Hb_Haemocue_observed = stan_data$Hb_Haemocue
  ind_haemocue = get_ind_haemocue(stan_data)
  Hb_Haemocue_mean_predicted = colMeans(Y_hat[,1,ind_haemocue])
  hist(Hb_Haemocue_observed-Hb_Haemocue_mean_predicted, breaks = 30,main='',
       xlab='Residuals (Haemocue haemoglobin)')
  boxplot((Hb_Haemocue_observed-Hb_Haemocue_mean_predicted) ~ stan_data$t_sim_hemocue_Hb)
  abline(h=0, lwd=2, col='red')
  
  ## Reticulocytes
  retic_observed = stan_data$Retic_data
  ind_retic = get_ind_retic(stan_data)
  retic_mean_predicted = colMeans(Y_hat[,2,ind_retic])
  hist(retic_observed-retic_mean_predicted, breaks = 30,
       xlab='Residuals (Reticulocyte count)',main='')
  boxplot((retic_observed-retic_mean_predicted) ~ stan_data$t_sim_retic)
  abline(h=0, lwd=2, col='red')
  
  pars = c(pars, 'theta_rand')
  
  my_thetas1 = rstan::extract(mod_fit_pop1, pars=pars)
  par(mfrow=c(2,2),las=1,family='serif',cex.lab=1.3, cex.axis=1.3)
  
  for(id in 1:23){
    ind_id = stan_data$ind_start_regimen[id]:stan_data$ind_end_regimen[id]
    Y_hat_hb = Y_hat[,1,ind_id]
    Y_hat_retic = Y_hat[,2,ind_id]
    
    ind_haemocue = stan_data$t_sim_hemocue_Hb[stan_data$ind_start_Hb_hemocue[id]:stan_data$ind_end_Hb_hemocue[id]]
    ind_hb_CBC = stan_data$t_sim_CBC_Hb[stan_data$ind_start_Hb_CBC[id]:stan_data$ind_end_Hb_CBC[id]]
    ind_retic = stan_data$t_sim_retic[stan_data$ind_start_retic[id]:stan_data$ind_end_retic[id]]
    
    plot(1:length(ind_id), colMeans(Y_hat_hb), type='l', lwd=3,
         panel.first=grid(), ylim =c(8.3,16.3), xlab='Time (days)',
         ylab='Haemoglobin (g/dL)',main=paste0('Subject ',id))
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
  effect = array(dim=c(length(my_thetas1$beta), length(xs)))
  for(i in 1:length(xs)){
    effect[,i] = emax_f(x = xs[i], 
                        logit_max_drug_effect = my_thetas1$logit_alpha,
                        emax_k = my_thetas1$h, 
                        ed50 = my_thetas1$beta)
  }
  effect_mean = colMeans(effect)
  plot(xs, 100*effect_mean,type='l',lwd=3, xlab='Primaquine mg/kg',
       ylab='Reduction in red cell lifespan (%)',panel.first=grid(),
       ylim = c(0,100))
  polygon(c(xs, rev(xs)),
          100*c(apply(effect,2, quantile,probs=0.1), rev(apply(effect,2,quantile,probs=0.9))),
          col=adjustcolor('grey',.3),border = NA)
  
  effect_id_list = list()
  for(id in 1:stan_data$N){
    logit_max_drug_effect = my_thetas1$logit_alpha+my_thetas1$theta_rand[,id,2]
    emax_k = my_thetas1$h
    ed50 = my_thetas1$beta*exp(my_thetas1$theta_rand[,id,3])
    
    
    effect_id = array(dim=c(length(my_thetas1$beta), length(xs)))
    for(i in 1:length(xs)){
      effect_id[,i]= emax_f(x = xs[i], 
                            logit_max_drug_effect = logit_max_drug_effect,
                            emax_k = emax_k, ed50 = ed50)
    }
    effect_id_list[[id]]=colMeans(effect_id)
    lines(xs, 100*colMeans(effect_id))
  }
  lines(xs, 100*colMeans(effect),lwd=3,col='red')
}



get_ind_cbc = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_CBC_Hb[stan_data$ind_start_Hb_CBC[i]:stan_data$ind_end_Hb_CBC[i]]])
  }
  return(ind)
}


get_ind_haemocue = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_hemocue_Hb[stan_data$ind_start_Hb_hemocue[i]:stan_data$ind_end_Hb_hemocue[i]]])
  }
  return(ind)
}

get_ind_retic = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_retic[stan_data$ind_start_retic[i]:stan_data$ind_end_retic[i]]])
  }
  return(ind)
}








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
