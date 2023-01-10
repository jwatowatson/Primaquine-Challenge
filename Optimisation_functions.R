## R code for the model
compute_effective_dose = function( drug_regimen,  nComp_sim,  K_weights,  mean_delay,  sigma_delay){
  
  lag=1;
  effective_dose = array(0, dim=length(drug_regimen));
  my_weights = array(dim = K_weights)
  for(t in 1:K_weights){
    my_weights[t] = pnorm(t, mean_delay, sigma_delay)-pnorm(t-1, mean_delay, sigma_delay);
  }
  for(t in 2:nComp_sim){
    if(t>K_weights){
      lag=lag+1;
      effective_dose[t] = sum(my_weights * drug_regimen[lag:t]);
    } else {      
      effective_dose[t] = sum(my_weights[(K_weights-t+1):K_weights] * drug_regimen[lag:t]);
    }
  }
  return (effective_dose);
}


dose_response = function(effective_dose, logit_max_drug_effect, emax_k, ed50){
  effect = gtools::inv.logit(logit_max_drug_effect) * (effective_dose^emax_k)/ (effective_dose^emax_k + ed50^emax_k);
  effect;
}


compute_rho = function(Hb, Hb_delta,  Hb_steady_state,
                       alpha_diff1,  alpha_delta1,
                       alpha_diff2,   alpha_delta2)
{
  rho_difference=0;
  rho_delta=0;
  rho=1; 
  
  if(Hb < Hb_steady_state) {
    # first increase in bone marrow production due to difference from steady state
    rho_difference= rho_difference+ alpha_diff1*(Hb_steady_state-Hb) + alpha_diff2*(Hb_steady_state-Hb)*(Hb_steady_state-Hb);
    # increase in bone marrow production due to fall in haemoglobin
    if(Hb_delta > 0) {
      rho_delta =rho_delta+ alpha_delta1*Hb_delta + alpha_delta2*Hb_delta*Hb_delta;
    }
  }
  # total increase in bone marrow production
  rho =rho+ rho_difference + rho_delta;
  return (rho);
}


# Compute the RBC age at which retics enter the circulation
# Hb: Hb: Hb level at time t
# Hb_steady_state: steady-state Hb
# Hb_transit_half: mid-point Hb at which the curve reaches half of the maximum transit time
compute_transit_time = function( Hb, Hb_steady_state, T_transit_steady_state, log_k)
{
  
  transit = T_transit_steady_state*exp(-(Hb_steady_state-Hb)*exp(log_k));
  if(transit>T_transit_steady_state) transit=T_transit_steady_state;
  transit;
}




# Count the number of retics in circulation
# transit: transit time
# reticulocytes: retics vector at time t
CountRetics = function(transit, reticulocytes){
  Total_retics = 0;
  T_retic = length(reticulocytes);
  i = T_retic;
  transit_int = ceiling(transit);
  retic_rem =  ceiling(transit) - transit;
  while(i > transit+1) {
    Total_retics = Total_retics + reticulocytes[i];
    i = i-1;
  }
  Total_retics = Total_retics + (retic_rem*reticulocytes[transit_int]);
  
  Total_retics;
}


# ************ deterministic simulation using the RBC dynamics model parameters **************
forwardsim = function(drug_regimen,            # dose given each day (mg/kg)
                      Hb_steady_state,           # Latent steady state value for the individual (g/dL)
                      alpha_diff1,               # maximum fold increase in RBC production
                      alpha_delta1,              # multiplicative factor on the fall in hb in increasing bone marrow production
                      alpha_diff2,               # maximum fold increase in RBC production
                      alpha_delta2,              # multiplicative factor on the fall in hb in increasing bone marrow production
                      logit_max_drug_effect,     # maximum drug effect
                      emax_k,          # slope effect in dose response function
                      ed50,                  # effective dose corresponding to half maximum drug effect (mg)
                      rbc_lifespan_steady_state, # steady-state RBC lifespan
                      log_k,                     # transit time parameter
                      nComp_sim,                 # Total number of compartments to forward simulate
                      T_nmblast,                 # normboblast lifespan - maturation time
                      T_retic,                   # reticulocyte lifespan - including marrow and circulating periods
                      T_RBC_max,                 # Maximum allowed value for RBC lifespan (arbitrary; days)
                      T_transit_steady_state,    # transit time at steady state
                      mean_delay,
                      sigma_delay,
                      K_weights,
                      sigma                      # sd for cumulative standard normal distribution in RBC survival function
){
  
  # set up variables we need for the simulation
  
  Hb_delta=0;
  
  effective_dose = array(dim=nComp_sim);       # The "effective dose" at each timepoint in the simulation
  Hb = array(nComp_sim);               # The haemoglobin at each point in the simulation
  retic_percent = array(nComp_sim);    # The retic count (%) at each point in the simulation
  temp_effect=0;
  
  # vectors to store red blood cell age distributions
  erythrocytes = temp_eryths = array(dim = T_RBC_max)
  normoblasts =  temp_normoblasts = array(dim = T_nmblast)
  reticulocytes = temp_retics = array(dim = T_retic);
  
  out_res = array(dim = c(2, nComp_sim))          # forming output of deterministic simulation
  
  # ***** Calculate initial values at steady state *****
  Hb[1] = Hb_steady_state;
  transit = compute_transit_time(Hb[1], Hb_steady_state, T_transit_steady_state, log_k);
  rho = compute_rho(Hb[1], Hb_delta, Hb_steady_state, alpha_diff1, alpha_delta1, alpha_diff2, alpha_delta2);
  
  # ***** Calculate the effective doses over time given the parameter mean_delay, sigma_delay
  effective_dose = compute_effective_dose(drug_regimen, nComp_sim, K_weights, mean_delay, sigma_delay);
  
  red_rbc_lifespan = rbc_lifespan_steady_state;
  lambda = 10000; # baseline production of normoblasts per unit time - arbitrary number
  
  erythrocytes[1] = lambda;
  for (i in 2:T_RBC_max){
    erythrocytes[i] = (1 - pnorm((i-red_rbc_lifespan)/sigma))*erythrocytes[i-1];
  }
  for(i in 1:T_nmblast){ normoblasts[i] = lambda; }
  for(i in 1:T_retic){ reticulocytes[i] = lambda; }
  
  Total_retics = CountRetics(transit, reticulocytes);
  Total_Eryths = sum(erythrocytes);
  C_0 = Total_retics + Total_Eryths;
  retic_percent[1] = 100 * Total_retics / C_0;
  
  #******************* Forward simulation *****************#
  for(t in 2:nComp_sim){
    
    # Compute the multiplication factor of basal normoblast production
    if(t>2) Hb_delta = Hb[t-2]-Hb[t-1];
    rho = compute_rho(Hb[t-1], Hb_delta, Hb_steady_state, alpha_diff1, alpha_delta1, alpha_diff2, alpha_delta2);
    
    # calculate the updated Hb dependent transit time for the reticulocytes
    transit = compute_transit_time(Hb[t-1], Hb_steady_state, T_transit_steady_state, log_k);
    
    # We move the RBCs from one compartment to the next
    temp_normoblasts = normoblasts;
    normoblasts[1] = rho*lambda;      # the number of new normoblasts made at time t
    for(i in 2:T_nmblast)             # move all the normoblasts along by one
    {
      normoblasts[i] = temp_normoblasts[i-1];
    }
    
    temp_retics = reticulocytes;
    reticulocytes[1] = temp_normoblasts[T_nmblast];
    for(i in 2:T_retic)
    {
      reticulocytes[i] = temp_retics[i-1];
    }
    
    # This part models the drug dependent killing as a shift in lifespan
    # Compute the life span at steady state for the current dose
    temp_effect = dose_response(effective_dose[t], logit_max_drug_effect, emax_k, ed50);
    red_rbc_lifespan = rbc_lifespan_steady_state*(1-temp_effect);
    
    
    # survival function of RBC lifespan
    temp_eryths = erythrocytes;
    erythrocytes[1] = temp_retics[T_retic];
    for (i in 2:T_RBC_max){
      erythrocytes[i] = (1 - pnorm((i-red_rbc_lifespan)/sigma))*temp_eryths[i-1];
    }
    
    # Count the number of retics and erythrocytes in circulation
    Total_retics = CountRetics(transit, reticulocytes);
    Total_Eryths = sum(erythrocytes);
    C_t = Total_retics + Total_Eryths;
    retic_percent[t] = 100 * Total_retics / C_t;
    Hb[t] = Hb_steady_state * C_t / C_0;
    
  }
  # store results
  out_res[1,] = Hb;
  out_res[2,] = retic_percent;
  
  out_res;
}

# test
drug_regimen = c(rep(10/60, 5), rep(15/60, 5), rep(30/60, 5), rep(0, 20))
pred=forwardsim(drug_regimen = drug_regimen,
                Hb_steady_state = 15, 
                alpha_diff1 = 0.01,
                alpha_delta1 = .5,
                alpha_diff2 = 0.003,
                alpha_delta2 = 0.15,
                logit_max_drug_effect = -0.45,
                emax_k = 1.7,
                ed50 = 0.2,
                rbc_lifespan_steady_state = 65,
                log_k = -2.3,
                nComp_sim = length(drug_regimen), 
                T_nmblast = 5,
                T_retic = 5, 
                T_RBC_max = 120, 
                T_transit_steady_state = 3.5, 
                mean_delay = 4, 
                sigma_delay = 1.9,
                K_weights = 10,
                sigma = 3)
# plot(pred[1,])
# plot(pred[2, ])
