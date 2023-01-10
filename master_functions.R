### make initialisation for stan ####
make_init_list = function(nchains){
  my_init = list()
  for(i in 1:nchains){
    my_init[[i]] = list(Hb_steady_state = rnorm(1, mean=14.5, sd=.1),
                        rbc_lifespan_steady_state = rnorm(1, mean = 75, sd = 5),
                        alpha_diff1 = 0.1,
                        alpha_diff2= 0.1,
                        alpha_delta1 = 0.1,
                        alpha_delta2 = 0.1,
                        logit_max_drug_effect = rnorm(1, mean = -0.6, sd = .1), # in logit scale
                        log_slope_effect = rnorm(1, mean = 1.3, sd = .11), # in log scale
                        log_ed50= rnorm(1, mean = log(15/60), sd = 1), # in log scale
                        log_k = rnorm(1, mean = -1.4, sd=.1), #in log scale
                        sigma_CBC = .5,
                        sigma_hemocue = .5,
                        sigma_retic_manual = .2,
                        sigma_Retic_data = .2) 
    names(my_init)[i] = paste('chain',i,sep='')
  }
  return(my_init)
}


###### stan data input #######

data_model = function(my_data, 
                      T_nmblast=5,
                      T_retic=5, 
                      T_RBC_max=140, 
                      T_transit_steady_state=3.5,
                      my_sigma = 3,
                      N_pred = NA,
                      drug_regimen_pred = NA){
  
  # check that properly sorted
  my_data = dplyr::arrange(my_data, ID, Study_Day)
  IDs = unique(my_data$ID)
  N = length(IDs)
  
  ind_start_regimen = ind_end_regimen = 
    ind_start_Hb_CBC = ind_end_Hb_CBC = 
    ind_start_Hb_hemocue = ind_end_Hb_hemocue = 
    ind_start_retic = ind_end_retic = 
    N_sim = 
    N_Hb_CBC = N_Hb_hemocue = N_retic = 
    array(dim=N)
  
  Hb_CBC = Hb_Haemocue = Retic_data = drug_regimen = 
    t_sim_CBC_Hb = t_sim_hemocue_Hb = t_sim_retic  = c()
  
  counter_Hb_CBC = counter_Hb_hemocue = counter_retic = counter_regimen = 1
  
  for(i in 1:N){
    
    ind_start_regimen[i] = counter_regimen
    
    ind_start_Hb_CBC[i] = counter_Hb_CBC 
    ind_start_Hb_hemocue[i] = counter_Hb_hemocue 
    ind_start_retic[i] = counter_retic

    
    id = IDs[i]
    data_id = dplyr::filter(my_data, ID==id) 
    
    N_sim[i] = max(data_id$Study_Day) - min(data_id$Study_Day) + 1 # max number of days to run forward simulation
    drug_regimen_id = array(0, dim = N_sim[i]) # initialise with all zeros
    
    timepoints = data_id$Study_Day - min(data_id$Study_Day) + 1
    drug_regimen_id[timepoints] = data_id$dosemgkg
    drug_regimen = c(drug_regimen, drug_regimen_id)
    
    ### --------------------------------
    ind_hemocue_Hb = !is.na(data_id$Haemocue_hb)
    ind_CBC_Hb = !is.na(data_id$CBC_hb)
    ind_retic = !is.na(data_id$Mean_Retic)

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

    ### --------------------------------
    Hb_haemocue_id= data_id$Haemocue_hb[ind_hemocue_Hb]
    Hb_CBC_id = data_id$CBC_hb[ind_CBC_Hb]
    Retic_id = data_id$Mean_Retic[ind_retic]

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
  
  N_CBC = length(Hb_CBC)
  N_hemo = length(Hb_Haemocue)
  N_retic = length(Retic_data)
  N_sim_tot = length(drug_regimen)
  
  data_stan = list(N=N, 
                   N_sim_tot=N_sim_tot, 
                   N_sim= N_sim, 
                   N_CBC = N_CBC, 
                   N_hemo = N_hemo, 
                   N_retic = N_retic, 
                   id = 1:N,
                   id2 = my_data$ID2[!duplicated(my_data$ID)],
                   drug_regimen = drug_regimen,
                   Hb_CBC = Hb_CBC,
                   Retic_data = round(Retic_data),
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
                   log_slope_effect_mean = 1, log_slope_effect_var = 1, #log scale
                   ed50_mean = 15/60, ed50_sd = .2,
                   prior_lifespan_SS_mean = 80, prior_lifespan_SS_var = 15,
                   Hb_steady_state_mean = 15, Hb_steady_state_var = 1,
                   log_k_mean = -1, log_k_var = 2, #log scale
                   max_drug_effect_mean = 0, max_drug_effect_var = 10, 
                   sigma_Hb_mean = .5, 
                   T_nmblast = T_nmblast, 
                   T_retic = T_retic, 
                   T_RBC_max = T_RBC_max,  
                   T_transit_steady_state = T_transit_steady_state,
                   sigma = my_sigma,
                   K_model=4,
                   N_pred=N_pred,
                   drug_regimen_pred=drug_regimen_pred)
  
  return(data_stan)
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
