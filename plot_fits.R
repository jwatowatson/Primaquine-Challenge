library(tidyverse)
library(gtools)
library(rstan)
source('master_functions.R')
source('Optimisation_functions.R')

pars_main = c('mean_delay','sigma_delay',
              'sigma_CBC','sigma_haemocue','sigma_retic',
              'logit_alpha','CBC_correction',
              'h','beta','log_k','Hb_star','T_E_star',
              'alpha_diff1','alpha_delta1' )

pars_RE = c('sigmasq_u','sigmasq_u_ic')

plot_model_fit = function(mod_fit, stan_data, pars) {
  
  traceplot(mod_fit, pars=pars)
  thetas_main = extract(mod_fit, pars=pars)
  hist(rweibull(10^4, shape = mean(thetas_main$mean_delay), 
                scale = thetas_main$sigma_delay))
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
  par(mfrow=c(3,2),las=1, family='serif')
  
  Y_hat = rstan::extract(mod_fit, pars = 'Y_hat')$Y_hat 
  ## CBC
  Hb_CBC_observed = stan_data$Hb_CBC
  ind_cbc = get_ind_cbc(stan_data)
  Hb_CBC_mean_predicted = colMeans(Y_hat[,1,ind_cbc]) + 
    mean(rstan::extract(mod_fit, pars='CBC_correction')$CBC_correction)
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
  
  data_pred_obs = tibble(
    day = c(stan_data$t_sim_CBC_Hb,
            stan_data$t_sim_hemocue_Hb),
    haemoglobin_obs = c(Hb_CBC_observed,
                        Hb_Haemocue_observed),
    haemoglobin_pred = c(Hb_CBC_mean_predicted,
                         Hb_Haemocue_mean_predicted)
  )
  med_vals = data_pred_obs %>% group_by(day) %>%
    summarise(Hb_obs = median(haemoglobin_obs),
              Hb_pred = median(haemoglobin_pred))
  plot(data_pred_obs$day,data_pred_obs$haemoglobin_obs,
       col=adjustcolor('grey',.4), panel.first=grid(), xlab='Day', ylab='Haemoglobin (g/dL)')
  lines(med_vals$day, med_vals$Hb_obs, lwd=3)
  lines(med_vals$day, med_vals$Hb_pred, lwd=3)
  
  
  pars = c(pars, 'theta_rand','theta_ic')
  
  my_thetas = rstan::extract(mod_fit, pars=pars)
  par(mfrow=c(4,4),las=1,family='serif',cex.lab=1.3, cex.axis=1.3)
  
  for(id in 1:stan_data$N_experiment){
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
  effect = array(dim=c(length(my_thetas$beta), length(xs)))
  for(i in 1:length(xs)){
    effect[,i] = emax_f(x = xs[i], 
                        logit_max_drug_effect = my_thetas$logit_alpha,
                        emax_k = my_thetas$h, 
                        ed50 = my_thetas$beta)
  }
  effect_mean = colMeans(effect)
  plot(xs, 100*effect_mean,type='l',lwd=3, xlab='Primaquine mg/kg',
       ylab='Reduction in red cell lifespan (%)',panel.first=grid(),
       ylim = c(0,100))
  polygon(c(xs, rev(xs)),
          100*c(apply(effect,2, quantile,probs=0.1), rev(apply(effect,2,quantile,probs=0.9))),
          col=adjustcolor('grey',.3),border = NA)
  
  effect_id_list = list()
  for(id in 1:stan_data$N_subject){
    logit_max_drug_effect = my_thetas$logit_alpha+my_thetas$theta_rand[,id,2]
    emax_k = my_thetas$h
    ed50 = my_thetas$beta*exp(my_thetas$theta_rand[,id,3])
    
    
    effect_id = array(dim=c(length(my_thetas$beta), length(xs)))
    for(i in 1:length(xs)){
      effect_id[,i]= emax_f(x = xs[i], 
                            logit_max_drug_effect = logit_max_drug_effect,
                            emax_k = emax_k, ed50 = ed50)
    }
    effect_id_list[[id]]=colMeans(effect_id)
    lines(xs, 100*colMeans(effect_id))
  }
  lines(xs, 100*colMeans(effect),lwd=3,col='red')
  
  
  
  
  
  
  drug_regimen = c(rep(15/60, 20), rep(0,15))
  eff_dose=compute_effective_dose(drug_regimen = drug_regimen,
                                  nComp_sim = length(drug_regimen),
                                  mean_delay = mean(thetas_main$mean_delay),
                                  sigma_delay = mean(thetas_main$sigma_delay))
  plot(1:length(drug_regimen), eff_dose, type='l')
}



get_ind_cbc = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N_experiment){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_CBC_Hb[stan_data$ind_start_Hb_CBC[i]:stan_data$ind_end_Hb_CBC[i]]])
  }
  return(ind)
}


get_ind_haemocue = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N_experiment){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_hemocue_Hb[stan_data$ind_start_Hb_hemocue[i]:stan_data$ind_end_Hb_hemocue[i]]])
  }
  return(ind)
}

get_ind_retic = function(stan_data){
  ind = c()
  for(i in 1:stan_data$N_experiment){
    ind = c(ind, (stan_data$ind_start_regimen[i]:stan_data$ind_end_regimen[i])[stan_data$t_sim_retic[stan_data$ind_start_retic[i]:stan_data$ind_end_retic[i]]])
  }
  return(ind)
}


get_rbc_lifespans = function(out, stan_data){
  thetas=rstan::extract(out, pars = c('T_E_star','theta_rand'))
  rbc_lifespans = array(dim=c(stan_data$N_experiment,dim(thetas$theta_rand)[1]))
  for(i in 1:stan_data$N_experiment){
    rbc_lifespans[i,] = thetas$theta_rand[,stan_data$id[i],8] + 
      thetas$T_E_star
    
  }
  return(rbc_lifespans)
}


get_predictions_LOO = function(out, stan_data){
  thetas = extract(out, pars = 'Y_pred')$Y_pred
  Y_pred_hb = data.frame(thetas[,1,])
  Y_pred_hb_rel = t(apply(Y_pred_hb, 1, function(x) 100*(1 - (x[1] - x)/x[1]) ))
  Y_pred_retic = data.frame(thetas[,2,])
  
  dat_true = data.frame(
    time = 1:length(stan_data$Y_true_haemocue),
    Y_true_haemocue=stan_data$Y_true_haemocue,
    Y_true_HbCBC=stan_data$Y_true_HbCBC,
    Y_true_retic=stan_data$Y_true_Retic)
  dat_true$baseline_hb = mean(c(dat_true$Y_true_haemocue[1], dat_true$Y_true_HbCBC[1]), na.rm = T)
  dat_true$Y_true_haemocue_rel = 100*(1 - (-dat_true$Y_true_haemocue + dat_true$baseline_hb)/dat_true$baseline_hb)
  dat_true$Y_true_CBC_rel = 100*(1 - (-dat_true$Y_true_HbCBC + dat_true$baseline_hb)/dat_true$baseline_hb)
  return(list(dat_true=dat_true, Y_pred_hb=Y_pred_hb, Y_pred_retic=Y_pred_retic, Y_pred_hb_rel=Y_pred_hb_rel))
}

load('Rout/stan_data_list.RData')

# 
pdf('model_fits_ep.pdf')
plot_model_fit(mod_fit = out_2,stan_data = dat_stan_list[[2]],
               pars = c(pars_main, pars_RE))
dev.off()



load('Rout/job_1.RData')
out_1 = out; check_rhat(out)

load('Rout/job_2.RData')
out_2 = out; check_rhat(out); rm(out)
rbcs = get_rbc_lifespans(out = out_2, stan_data = dat_stan_list[[2]])
plot(rowMeans(rbcs))
plot(out_2, pars='dose_weights')
traceplot(out_1, pars=pars_main)
traceplot(out_2, pars=pars_main)


pdf('LOO_predictions.pdf', width = 10, height = 10)
par(las=1, mfrow=c(4,4), family='serif')

for(kk in 1:dat_stan_list[[1]]$N_experiment){
  load(paste0('Rout/job_',kk+2,'.RData'))
  # rstan::traceplot(out,pars=pars_main)
  preds_LOO=get_predictions_LOO(out = out, stan_data = dat_stan_list[[kk+2]])
  
  ## Haemolgobin
  plot(preds_LOO$dat_true$time, preds_LOO$dat_true$Y_true_haemocue_rel,pch=15,
       ylab='Change from baseline Hb (%)', xlab = 'days', 
       panel.first=grid(), ylim = c(60,110))
  points(preds_LOO$dat_true$Y_true_CBC_rel, pch=16)
  title(unique(PQdat$ID2)[kk])
  
  lines(preds_LOO$dat_true$time, colMeans(preds_LOO$Y_pred_hb_rel),lwd=3)
  for(j in 1:100){
    k=sample(size = 1, x = nrow(preds_LOO$Y_pred_hb))
    lines(preds_LOO$dat_true$time, preds_LOO$Y_pred_hb_rel[k,], col=adjustcolor('grey',.3))
  }
  points(preds_LOO$dat_true$Y_true_haemocue_rel, pch=15,col='red')
  points(preds_LOO$dat_true$Y_true_CBC_rel, pch=16,col='red')
  
  ## Reticulocytes
  plot(preds_LOO$dat_true$time, preds_LOO$dat_true$Y_true_retic,pch=15,
       ylab='Reticulocytes (%)', xlab = 'days', 
       panel.first=grid(), ylim = c(0,15))
  title(unique(PQdat$ID2)[kk])
  
  lines(preds_LOO$dat_true$time, colMeans(preds_LOO$Y_pred_retic),lwd=3)
  for(j in 1:100){
    k=sample(size = 1, x = nrow(preds_LOO$Y_pred_hb))
    lines(preds_LOO$dat_true$time, preds_LOO$Y_pred_retic[k,], col=adjustcolor('grey',.3))
  }
  points(preds_LOO$dat_true$Y_true_retic, pch=15,col='red')
}

dev.off()

job_files = list.files('Rout', pattern = 'job', full.names = T)
preds_list = list()
for(kk in 3:length(job_files)){
  load(paste0('Rout/job_',kk,'.RData'))
  rhat_vals = summary(out)$summary[,'Rhat']
  if(max(rhat_vals)>1.2) {writeLines(sprintf('Max Rhat is %s', max(rhat_vals)))
  } else { 
    writeLines('Convergence OK')
  }
  thetas = extract(out)
}

plot(out_3, pars='dose_weights')

thetas = extract(out, pars='Y_hat')$Y_hat[,3,]
plot(colMeans(thetas))

par(mfrow=c(2,2))
for(id in dat_stan_list[[3]]$id){
  ind = dat_stan_list[[3]]$ind_start_regimen[id]:dat_stan_list[[3]]$ind_end_regimen[id]
  plot(dat_stan_list[[3]]$drug_regimen[ind],type='l',lwd=3, ylim = c(0,1))
  lines(colMeans(thetas)[ind])
}


