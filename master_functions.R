### make initialisation for stan ####
make_init_list = function(nchains){
  my_init = list()
  for(i in 1:nchains){
    my_init[[i]] = list(Hb_star = rnorm(1, mean=15, sd=.5),
                        T_E_star = rnorm(1, mean = 90, sd = 10),
                        alpha_diff1 = rexp(1, rate = 40),
                        alpha_diff2 = rexp(1, rate = 40),
                        alpha_delta1 = rexp(1, rate = 10),
                        alpha_delta2 = rexp(1, rate = 10),
                        diff_alpha = 0.05,
                        delta_alpha=0.05,
                        logit_alpha = rnorm(1, mean = -0.5, sd = .25), # on logit scale
                        beta = runif(1, min = 0.1, max = 0.4),
                        h = rexp(1),
                        log_k = rnorm(1, mean = -2, sd=.5), #in log scale
                        mean_delay = runif(1, min = 1, max = 10),
                        sigma_delay = 1+rexp(1, rate = 2),
                        sigma_CBC = rexp(1),
                        sigma_manual = rexp(1),
                        cumdose_alpha=c(1,3),
                        logit_MAX_EFFECT = -2.5,
                        logit_G6PD_delta_day= -3,
                        logit_G6PD_threshold= -5,
                        CBC_correction = rnorm(1)) 
    names(my_init)[i] = paste('chain',i,sep='')
  }
  return(my_init)
}


###### stan data input #######
make_stan_dataset = function(my_data, 
                             ID_experiment = 'ID2',
                             ID_subject = 'ID',
                             T_nmblast=5,
                             T_retic=5, 
                             T_RBC_max=140, 
                             T_transit_steady_state=3.5,
                             my_sigma = 3,
                             N_pred = NA,
                             K_weights = 10,
                             drug_regimen_pred = NA,
                             data_pred = NA){
  
  
  # check that properly sorted
  key_cols = c('Study_Day','dosemgkg','Haemocue_hb','CBC_hb','Mean_retic')
  if(!all(key_cols %in% colnames(my_data))) {
    stop(sprintf('missing key column %s', key_cols[!key_cols %in% colnames(my_data)]))
  }
  
  IDs_experiment = unique(my_data[, ID_experiment,drop=T])
  IDs_subject = unique(my_data[, ID_subject,drop=T])
  N_subject = length(IDs_subject)
  N_experiment = length(IDs_experiment)
  
  my_data$ID_experiment = as.numeric(as.factor(my_data[, ID_experiment, drop=T]))
  my_data$ID_experiment_char = my_data[, ID_experiment, drop=T]
  
  my_data$ID_subject = as.numeric(as.factor(my_data[, ID_subject, drop=T]))
  my_data$ID_subject_char = my_data[, ID_subject, drop=T]
  
  ind_start_regimen = ind_end_regimen = 
    ind_start_Hb_CBC = ind_end_Hb_CBC = 
    ind_start_Hb_hemocue = ind_end_Hb_hemocue = 
    ind_start_retic = ind_end_retic = 
    N_sim = 
    N_Hb_CBC = N_Hb_hemocue = N_retic = 
    array(dim=N_experiment)
  
  Hb_CBC = Hb_Haemocue = Retic_data = drug_regimen = 
    t_sim_CBC_Hb = t_sim_hemocue_Hb = t_sim_retic  = 
    t_sim_CBC_Hb_study = t_sim_hemocue_Hb_study = t_sim_retic_study = c()
  
  counter_Hb_CBC = counter_Hb_hemocue = counter_retic = counter_regimen = 1
  
  ID_random_effect = array(dim = N_experiment)
  ID_repeat_random_effect = array(0, dim = N_experiment); k_repeat=0
  
  my_data$study_numeric = as.numeric(as.factor(my_data$study))
  for(i in unique(my_data$ID_experiment)){
    
    ind_start_regimen[i] = counter_regimen
    
    ind_start_Hb_CBC[i] = counter_Hb_CBC 
    ind_start_Hb_hemocue[i] = counter_Hb_hemocue 
    ind_start_retic[i] = counter_retic
    
    data_exp = my_data %>% filter(ID_experiment==i) 
    ID_random_effect[i] = unique(data_exp$ID_subject)
    if(duplicated(ID_random_effect)[i]) {
      k_repeat = k_repeat+1
      ID_repeat_random_effect[i] = k_repeat
    }
    
    N_sim[i] = max(data_exp$Study_Day) - min(data_exp$Study_Day) + 1 # max number of days to run forward simulation
    drug_regimen_id = array(0, dim = N_sim[i])     # initialise with all zeros
    
    timepoints = data_exp$Study_Day - min(data_exp$Study_Day) + 1 ## observed days
    drug_regimen_id[timepoints] = data_exp$dosemgkg
    drug_regimen = c(drug_regimen, drug_regimen_id)
    
    ### --------------------------------
    ind_hemocue_Hb = !is.na(data_exp$Haemocue_hb)
    ind_CBC_Hb = !is.na(data_exp$CBC_hb)
    ind_retic = !is.na(data_exp$Mean_retic)
    
    N_Hb_hemocue[i] = sum(ind_hemocue_Hb)
    N_Hb_CBC[i] = sum(ind_CBC_Hb)
    N_retic[i] = sum(ind_retic)
    
    ### --------------------------------
    t_sim_hemocue_Hb_id = timepoints[ind_hemocue_Hb]
    t_sim_CBC_Hb_id = timepoints[ind_CBC_Hb]
    t_sim_retic_id = timepoints[ind_retic]
    
    t_sim_hemocue_Hb = c(t_sim_hemocue_Hb, t_sim_hemocue_Hb_id)
    t_sim_CBC_Hb = c(t_sim_CBC_Hb, t_sim_CBC_Hb_id)
    t_sim_retic = c(t_sim_retic, t_sim_retic_id)
    
    t_sim_hemocue_Hb_study = c(t_sim_hemocue_Hb_study,data_exp$study_numeric[ind_hemocue_Hb])
    t_sim_CBC_Hb_study = c(t_sim_CBC_Hb_study, data_exp$study_numeric[ind_CBC_Hb])
    t_sim_retic_study = c(t_sim_retic_study, data_exp$study_numeric[ind_retic])
    
    ### --------------------------------
    Hb_haemocue_id= data_exp$Haemocue_hb[ind_hemocue_Hb]
    Hb_CBC_id = data_exp$CBC_hb[ind_CBC_Hb]
    Retic_id = data_exp$Mean_retic[ind_retic]
    
    Hb_CBC = c(Hb_CBC, Hb_CBC_id)
    Hb_Haemocue = c(Hb_Haemocue, Hb_haemocue_id)
    Retic_data = c(Retic_data, Retic_id)
    
    ### --------------------------------
    ind_end_regimen[i] = counter_regimen+N_sim[i]-1
    counter_regimen=counter_regimen+N_sim[i]
    
    ind_end_Hb_CBC[i] = counter_Hb_CBC + N_Hb_CBC[i]-1
    ind_end_Hb_hemocue[i] = counter_Hb_hemocue + N_Hb_hemocue[i]-1
    ind_end_retic[i] = counter_retic + N_retic[i]-1
    
    counter_Hb_CBC = counter_Hb_CBC + N_Hb_CBC[i]
    counter_Hb_hemocue = counter_Hb_hemocue + N_Hb_hemocue[i]
    counter_retic = counter_retic + N_retic[i]
  }
  
  if(is.na(N_pred) | is.na(drug_regimen_pred)){
    drug_regimen_pred=rep(0.25, 21)
    N_pred = length(drug_regimen_pred)
    Y_true_haemocue=NA
    Y_true_HbCBC=NA
    Y_true_Retic=NA
  } 
  if (!all(is.na(data_pred))) {
    t_pred = data_pred$Study_Day
    N_pred = max(data_pred$Study_Day)-min(data_pred$Study_Day)+1
    drug_regimen_pred = rep(0, N_pred)
    drug_regimen_pred[data_pred$Study_Day-min(data_pred$Study_Day)+1]=data_pred$dosemgkg
    Y_true_haemocue = rep(NA, N_pred)
    Y_true_haemocue[data_pred$Study_Day-min(data_pred$Study_Day)+1]=data_pred$Haemocue_hb
    Y_true_HbCBC = rep(NA, N_pred)
    Y_true_HbCBC[data_pred$Study_Day-min(data_pred$Study_Day)+1]=data_pred$CBC_hb
    Y_true_Retic = rep(NA, N_pred)
    Y_true_Retic[data_pred$Study_Day-min(data_pred$Study_Day)+1]=data_pred$Manual_retic
  }
  
  N_CBC = length(Hb_CBC)
  N_hemo = length(Hb_Haemocue)
  N_retic = length(Retic_data)
  N_sim_tot = length(drug_regimen)
  
  data_stan = list(N_experiment=N_experiment, 
                   N_subject = N_subject,
                   N_sim_tot=N_sim_tot, 
                   N_sim= N_sim, 
                   N_CBC = N_CBC, 
                   N_hemo = N_hemo, 
                   N_retic = N_retic,
                   N_repeat = k_repeat,
                   repeat_oc = ID_repeat_random_effect,
                   id = ID_random_effect,
                   drug_regimen = drug_regimen,
                   Hb_CBC = Hb_CBC,
                   Retic_data = Retic_data,
                   Hb_Haemocue = Hb_Haemocue, 
                   ind_start_regimen = ind_start_regimen,
                   ind_end_regimen = ind_end_regimen,
                   ind_start_Hb_CBC = ind_start_Hb_CBC,
                   ind_end_Hb_CBC = ind_end_Hb_CBC, 
                   ind_start_Hb_hemocue = ind_start_Hb_hemocue, 
                   ind_end_Hb_hemocue = ind_end_Hb_hemocue, 
                   ind_start_retic = ind_start_retic, 
                   ind_end_retic = ind_end_retic, 
                   t_sim_hemocue_Hb = t_sim_hemocue_Hb, 
                   t_sim_CBC_Hb = t_sim_CBC_Hb, 
                   t_sim_retic = t_sim_retic, 
                   t_sim_hemocue_Hb_study=t_sim_hemocue_Hb_study,
                   t_sim_CBC_Hb_study=t_sim_CBC_Hb_study,
                   t_sim_retic_study=t_sim_retic_study,
                   log_slope_effect_mean = 1, log_slope_effect_var = 1, #log scale
                   beta_mean = 15/60, 
                   beta_sigma = .2,
                   T_E_star_mean = 90, 
                   T_E_star_sigma = 20,
                   Hb_star_mean = 15,
                   Hb_star_sigma = 1,
                   log_k_mean = -1, 
                   log_k_sigma = 2, #log scale
                   alpha_mean = 0, 
                   alpha_sigma = 2, 
                   sigma_Hb_mean = 2, 
                   T_nmblast = T_nmblast, 
                   T_retic = T_retic, 
                   T_RBC_max = T_RBC_max,  
                   T_transit_steady_state = T_transit_steady_state,
                   sigma = my_sigma,
                   N_pred=N_pred,
                   drug_regimen_pred=drug_regimen_pred,
                   K_weights = K_weights,
                   Y_true_haemocue=Y_true_haemocue,
                   Y_true_HbCBC=Y_true_HbCBC,
                   Y_true_Retic=Y_true_Retic)
  
  return(data_stan)
}


# check Rhat

check_rhat = function(out){
  rhat_vals = summary(out)$summary[,'Rhat']
  return(max(rhat_vals))
}



###### FUNCTION: Extract and format Stan output ######
extract_draws = function(parExtract, fit) {
  fit %>%
    rstan::extract(., pars = parExtract) %>%
    as.data.frame(.) %>%
    gather(key = "par", value = "value") %>% .$value
}


## EMAX function ####
emax_f = function(x, logit_max_drug_effect, emax_k,ed50 ){
  out = (gtools::inv.logit(logit_max_drug_effect) * x^emax_k )/ (x^emax_k + ed50^emax_k);
  return(out)
}
