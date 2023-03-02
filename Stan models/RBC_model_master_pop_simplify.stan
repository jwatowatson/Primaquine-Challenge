//
//###############################################################################
//# Authors: James Watson and Parinaz Mehdipour                                 #
//# Date: March 2022                                                            #
//# Name: RBC model_master_pop.stan                                             #
//# Description: Stan code for red blood cell model                             #
//###############################################################################
//# This code includes two parts including the red blood cell dynamics functions
//# and the hierarchial model fitting

//Change the alpha_diff and alpha_delta parameters to the linear scale (not log) – hence constrain to be >0
//Make the random effects for alpha_diff/delta to be on log scale (so individual parameter is eg alpha_diff * exp(theta_id) )
//Put normal priors on the alpha_diff/delta random effects variance terms (so you can’t get massive variance estimates)


functions {
  // This function computes the effective dose given a dosing schedule vector
  // drug_regimen: a vector of daily primaquine / `equivalent' doses in mg
  // (in practice only dosed once a day, this can be changed to hours, weeks, ...)
  // t: current time point
  vector compute_effective_dose(vector drug_regimen, int nComp_sim, int K_weights, real mean_delay, real sigma_delay){
    vector[nComp_sim] effective_dose;
    vector[K_weights] weights;
    int lag=1;
    effective_dose[1]=0;
    
    // make the weights schema
    for(t in 1:K_weights){
      weights[t] = normal_cdf(t, mean_delay, sigma_delay)-normal_cdf(t-1, mean_delay, sigma_delay);
    }
    weights = weights/sum(weights);// normalise so they sum to 1
    
    for(t in 2:nComp_sim){
      if(t>K_weights){
        lag+=1;
        effective_dose[t] = sum(weights .* drug_regimen[lag:t]);
      } else {      
        effective_dose[t] = sum(weights[(K_weights-t+1):K_weights] .* drug_regimen[lag:t]);
      }
    }
    return effective_dose;
  }
  
  
  // This function computes the reduction in RBC lifespan
  // effective_dose: effective dose by current time t
  // logit_alpha: maximum % decrease in RBC lifespan - logit scale
  // h: slope of the curve in log scale
  // beta: effective dose value when maximum drug effect is halved.
  real dose_response(real effective_dose, real logit_alpha, real h, real beta){
    real effect;
    effect = inv_logit(logit_alpha) * pow(effective_dose,h)/ (pow(effective_dose,h) + pow(beta,h));
    return effect;
  }
  
  // Compute the fold change from baseline in the number of new normoblasts produced given the Hb level
  // Hb: Hb level at time t
  // Hb_star: steady-state Hb
  // rho_max: maximum fold increase in production relative to steady state
  // Hb_rho_half: mid-point Hb when maximum production is halved
  real compute_rho(real Hb, real Hb_delta, real Hb_star,
  real alpha_diff1, real alpha_delta1,
  real alpha_diff2,  real alpha_delta2)
  {
    real rho_difference=0;
    real rho_delta=0;
    real rho=1; // rho = 1 is baseline production
    
    if(Hb < Hb_star) {
      // first increase in bone marrow production due to difference from steady state
      rho_difference += alpha_diff1*(Hb_star-Hb) + alpha_diff2*(Hb_star-Hb)*(Hb_star-Hb);
      // increase in bone marrow production due to fall in haemoglobin
      if(Hb_delta > 0) {
        rho_delta += alpha_delta1*Hb_delta + alpha_delta2*Hb_delta*Hb_delta;
      }
    }
    // total increase in bone marrow production
    rho += rho_difference + rho_delta;
    return rho;
  }
  
  
  // Compute the RBC age at which retics enter the circulation
  // Hb: Hb: Hb level at time t
  // Hb_star: steady-state Hb
  // Hb_transit_half: mid-point Hb at which the curve reaches half of the maximum transit time
  real compute_transit_time(real Hb, real Hb_star, real T_transit_steady_state, real log_k)
  {
    real transit;
    transit = T_transit_steady_state*exp(-(Hb_star-Hb)*exp(log_k));
    if(transit>T_transit_steady_state) transit=T_transit_steady_state;
    return transit;
  }
  
  
  //function for changing real to integer (for transit time)
  int real2int(real x, int min_val, int max_val){
    int out;
    real y = round(x);
    for (n in min_val:max_val) {
      if (n >= y){
        out = n;
        return out;
      }
    }
    return out;
  }
  
  
  // Count the number of retics in circulation
  // transit: transit time
  // reticulocytes: retics vector at time t
  real CountRetics(real transit, vector reticulocytes){
    real Total_retics = 0;
    int T_retic = rows(reticulocytes);
    int i = T_retic;
    real transit_idx = ceil(transit);
    int transit_int = real2int(transit_idx, 0, 6);
    real retic_rem =  ceil(transit) - transit;
    while(i > transit+1) {
      Total_retics += reticulocytes[i];
      i = i-1;
    }
    Total_retics += (retic_rem*reticulocytes[transit_int]);
    
    return Total_retics;
  }
  
  
  // ************ deterministic simulation using the RBC dynamics model parameters **************
  
  matrix forwardsim(vector drug_regimen,            // dose given each day (mg/kg)
  real Hb_star,           // Latent steady state value for the individual (g/dL)
  real alpha_diff1,               // maximum fold increase in RBC production
  real alpha_delta1,              // multiplicative factor on the fall in hb in increasing bone marrow production
  real alpha_diff2,               // maximum fold increase in RBC production
  real alpha_delta2,              // multiplicative factor on the fall in hb in increasing bone marrow production
  real logit_alpha,     // maximum drug effect
  real h,          // slope effect in dose response function
  real beta,                  // effective dose corresponding to half maximum drug effect (mg)
  real T_E_star, // steady-state RBC lifespan
  real log_k,                     // transit time parameter
  int  nComp_sim,                 // Total number of compartments to forward simulate
  int  T_nmblast,                 // normboblast lifespan - maturation time
  int  T_retic,                   // reticulocyte lifespan - including marrow and circulating periods
  int  T_RBC_max,                 // Maximum allowed value for RBC lifespan (arbitrary; days)
  real T_transit_steady_state,    // transit time at steady state
  real mean_delay,
  real sigma_delay,
  int K_weights,
  real sigma                      // sd for cumulative standard normal distribution in RBC survival function
  ){
    
    // set up variables we need for the simulation
    real transit;
    real red_rbc_lifespan;
    real rho;
    real Hb_delta=0;
    real lambda;                               // baseline production of normoblasts per time
    real Total_Eryths;                         // Total circulating erythrocytes
    real Total_retics;                         // Total circulating reticulocytes
    
    vector[nComp_sim] effective_dose;       // The "effective dose" at each timepoint in the simulation
    row_vector[nComp_sim] Hb;               // The haemoglobin at each point in the simulation
    row_vector[nComp_sim] retic_percent;    // The retic count (%) at each point in the simulation
    real temp_effect=0;
    real C_t;                               // Counter of cells
    real C_0;                               // Counter of cells - at baseline steady state
    
    // vectors to store red blood cell age distributions
    vector[T_RBC_max] erythrocytes;         // Distribution of erythrocytes in circulation
    vector[T_RBC_max] temp_eryths;
    vector[T_nmblast] normoblasts;          // Distribution of normoblasts
    vector[T_nmblast] temp_normoblasts;
    vector[T_retic] reticulocytes;          // Distribution of reticulocytes
    vector[T_retic] temp_retics;
    
    matrix[2, nComp_sim] out_res;           // forming output of deterministic simulation
    
    // ***** Calculate initial values at steady state *****
    Hb[1] = Hb_star;
    transit = compute_transit_time(Hb[1], Hb_star, T_transit_steady_state, log_k);
    rho = compute_rho(Hb[1], Hb_delta, Hb_star, alpha_diff1, alpha_delta1, alpha_diff2, alpha_delta2);
    
    // ***** Calculate the effective doses over time given the parameter mean_delay, sigma_delay
    effective_dose = compute_effective_dose(drug_regimen, nComp_sim, K_weights, mean_delay, sigma_delay);
    
    red_rbc_lifespan = T_E_star;
    lambda = 10000; // baseline production of normoblasts per unit time - arbitrary number
    
    erythrocytes[1] = lambda;
    for (i in 2:T_RBC_max){
      erythrocytes[i] = (1 - Phi_approx((i-red_rbc_lifespan)/sigma))*erythrocytes[i-1];
    }
    for(i in 1:T_nmblast){ normoblasts[i] = lambda; }
    for(i in 1:T_retic){ reticulocytes[i] = lambda; }
    
    Total_retics = CountRetics(transit, reticulocytes);
    Total_Eryths = sum(erythrocytes);
    C_0 = Total_retics + Total_Eryths;
    retic_percent[1] = 100 * Total_retics / C_0;
    
    //******************* Forward simulation *****************//
    for(t in 2:nComp_sim){
      
      // Compute the multiplication factor of basal normoblast production
      if(t>2) Hb_delta = Hb[t-2]-Hb[t-1];
      rho = compute_rho(Hb[t-1], Hb_delta, Hb_star, alpha_diff1, alpha_delta1, alpha_diff2, alpha_delta2);
      
      // calculate the updated Hb dependent transit time for the reticulocytes
      transit = compute_transit_time(Hb[t-1], Hb_star, T_transit_steady_state, log_k);
      
      // We move the RBCs from one compartment to the next
      temp_normoblasts = normoblasts;
      normoblasts[1] = rho*lambda;      // the number of new normoblasts made at time t
      for(i in 2:T_nmblast)             // move all the normoblasts along by one
      {
        normoblasts[i] = temp_normoblasts[i-1];
      }
      
      temp_retics = reticulocytes;
      reticulocytes[1] = temp_normoblasts[T_nmblast];
      for(i in 2:T_retic)
      {
        reticulocytes[i] = temp_retics[i-1];
      }
      
      // This part models the drug dependent killing as a shift in lifespan
      // Compute the life span at steady state for the current dose
      temp_effect = dose_response(effective_dose[t], logit_alpha, h, beta);
      red_rbc_lifespan = T_E_star*(1-temp_effect);
      
      // survival function of RBC lifespan
      temp_eryths = erythrocytes;
      erythrocytes[1] = temp_retics[T_retic];
      for (i in 2:T_RBC_max){
        erythrocytes[i] = (1 - Phi_approx((i-red_rbc_lifespan)/sigma))*temp_eryths[i-1];
      }
      
      // Count the number of retics and erythrocytes in circulation
      Total_retics = CountRetics(transit, reticulocytes);
      Total_Eryths = sum(erythrocytes);
      C_t = Total_retics + Total_Eryths;
      retic_percent[t] = 100 * Total_retics / C_t;
      Hb[t] = Hb_star * C_t / C_0;
      
    }
    // store results
    out_res[1] = Hb;
    out_res[2] = retic_percent;
    
    return out_res;
  }
}

//////////////////////////////////////////////////////////
/////////////////////////Fitting Model////////////////////
/////////////////////////////////////////////////////////
data {
  
  // fixed parameters for the model
  int T_nmblast;                 // normboblast lifespan - maturation time
  int T_retic;                   // reticulocyte lifespan - including marrow and circulating periods
  int T_RBC_max;                 // Maximum allowed value for RBC lifespan (arbitrary; days)
  real T_transit_steady_state;   // 
  real sigma;                    // sd for cumulative standard normal distribution in RBC survival function
  int K_weights;

  // Defining data to be fitted
  int<lower=1> N;          // number of subjects
  int<lower=0> N_sim_tot;  // total number of days to simulate across all patients
  int<lower=0> N_CBC;      // total number of Hb - CBC observations
  int<lower=0> N_hemo;     // total number of Hb - Haemocue observations
  int<lower=0> N_retic;    // total number of retic (average of CBC and manual readings)

  // observed data and indexing
  // drug regimen to use and indexing
  int<lower=1> N_sim[N];              // length of simulation for each individual
  int<lower=1> ind_start_regimen[N];  // starting indices for each individual
  int<lower=1> ind_end_regimen[N];    // ending indices for each individual
  vector[N_sim_tot] drug_regimen;     // doses array (all individuals)
  int<lower=1,upper=N> id[N];         // vector for subject ID (used if repeated experiments)
  
  // - Hb CBC
  int<lower=1> ind_start_Hb_CBC[N];
  int<lower=1> ind_end_Hb_CBC[N];
  int<lower=1> t_sim_CBC_Hb[N_CBC];   // days of observation
  real<lower=0> Hb_CBC[N_CBC];
  
  // - Hb manual read
  int<lower=1> ind_start_Hb_hemocue[N];
  int<lower=1> ind_end_Hb_hemocue[N];
  int<lower=1> t_sim_hemocue_Hb[N_hemo];   // days with observed data
  real<lower=0> Hb_Haemocue[N_hemo];
  
  // - Reticulocyte data
  int<lower=1> ind_start_retic[N];
  int<lower=1> ind_end_retic[N];
  int<lower=1> t_sim_retic[N_retic];   // days of observation
  real<lower=0> Retic_data[N_retic];
  
  // Prior parameters
  //mean
  real beta_mean;
  real T_E_star_mean;
  real<lower=0> Hb_star_mean;
  real log_k_mean;
  real alpha_mean;
  
  //variation
  real<lower=0> log_slope_effect_var;
  real<lower=0> beta_sigma;
  real<lower=0> T_E_star_sigma;
  real<lower=0> Hb_star_sigma;
  real<lower=0> log_k_sigma;
  real<lower=0> alpha_sigma;
  
  //residual
  real<lower=0> sigma_Hb_mean;
  
  // single drug regimen to predict
  int<lower=1> N_pred;
  vector[N_pred] drug_regimen_pred;
}

transformed data{
  int K_rand_effects = 4;
  vector[K_rand_effects] my_zeros;
  for(i in 1:K_rand_effects) my_zeros[i] = 0;
}

parameters {
  
  // error terms
  real<lower=0> sigma_CBC;
  real<lower=0> sigma_haemocue;
  real<lower=0> sigma_retic;
  real CBC_correction; // systematic difference between CBC haemoglobin and Haemocue haemoglobin
  
  // parameters governing the dose-response curve
  real logit_alpha; // maximal effect on logit scale
  real<lower=0> beta; // half maximal effect
  real<lower=0> h; // hill coefficient
  
  // parameters governing the delay in effect
  real<lower=1,upper=10> mean_delay;
  real<lower=0> sigma_delay; 
  
  // retic transit function
  real log_k;
  
  real<lower=0> Hb_star; // steady state haemoglobin
  real<lower=0> T_E_star;// steady state red cell lifespan
  
  // parameters governing the production of new cells in bone marrow
  real<lower=0> alpha_diff1;
  real<lower=0> alpha_delta1;
  real<lower=0> alpha_diff2;
  real<lower=0> alpha_delta2;
  
  // random effects
  vector[K_rand_effects] theta_rand[N]; // individual random effects vector
  cholesky_factor_corr[K_rand_effects] L_Omega; // correlation matrix
  vector<lower=0>[K_rand_effects] sigmasq_u; // variance of random effects
}

transformed parameters {
  matrix[2,N_sim_tot] Y_hat;
  
  for(j in 1:N){
    
    // compute output from forward_sim
    Y_hat[,ind_start_regimen[j]:ind_end_regimen[j]] = forwardsim(drug_regimen[ind_start_regimen[j]:ind_end_regimen[j]],
    Hb_star+theta_rand[id[j]][1], // steady state haemoglobin
    alpha_diff1, // parameters on bone marrow response (polynomial)
    alpha_delta1,
    alpha_diff2,
    alpha_delta2,
    logit_alpha+theta_rand[id[j]][2],// max effect on lifespan
    h, // slope of effect on lifespan
    beta*exp(theta_rand[id[j]][3]),// dose giving half max effect on lifespan
    T_E_star+theta_rand[id[j]][4],// steady state lifespan
    log_k,// retic release parameter
    N_sim[j],
    T_nmblast, T_retic, T_RBC_max, T_transit_steady_state, 
    mean_delay,
    sigma_delay,
    K_weights,
    sigma);
  }
}

model{
  // Prior
  // error terms
  sigma_CBC ~ exponential(1);
  sigma_haemocue ~ exponential(1);
  sigma_retic ~ normal(0.18, .1);
  CBC_correction ~ normal(0,1); 
  
  // random effects
  sigmasq_u[1] ~ exponential(2);   // prior on variance: steady state hb
  sigmasq_u[2] ~ exponential(2);   // prior: logit_alpha
  sigmasq_u[3] ~ exponential(2);   // prior: beta
  sigmasq_u[4] ~ exponential(1);   // prior: rbc lifespan

  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:N) theta_rand[i] ~ multi_normal_cholesky(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
  
  // steady state parameters
  T_E_star ~ normal(T_E_star_mean,T_E_star_sigma);
  Hb_star ~ normal(Hb_star_mean,Hb_star_sigma);
  
  // delay in dose effect
  mean_delay ~ normal(5,5) T[0,10];
  sigma_delay ~ exponential(1);
  
  // parameters governing the dose-response curve
  h ~ exponential(1);
  logit_alpha ~ normal(alpha_mean, alpha_sigma);
  beta ~ normal(beta_mean, beta_sigma) T[0,];
  
  // parameters governing the production of new cells in bone marrow
  alpha_delta1 ~ exponential(1);
  alpha_diff1 ~  exponential(1);
  alpha_delta2 ~ exponential(1);
  alpha_diff2 ~  exponential(1);
  
  // retic transit function
  log_k ~ normal(log_k_mean,log_k_sigma);
  
  // Likelihood
  for(j in 1:N){
    // calculate likelihood
    Hb_CBC[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]] ~ student_t(10,Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_CBC_Hb[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]]]+CBC_correction, sigma_CBC);
    Hb_Haemocue[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]] ~ student_t(10,Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_hemocue_Hb[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]]], sigma_haemocue);
    log(Retic_data[ind_start_retic[j]:ind_end_retic[j]]) ~ student_t(10, log(Y_hat[2][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_retic[ind_start_retic[j]:ind_end_retic[j]]]), sigma_retic);
  }
  
}


generated quantities {
  
  matrix[2,N_pred] Y_pred;
  {
    vector[K_rand_effects] theta_rand_pred; // individual random effects vector
  
    theta_rand_pred = multi_normal_cholesky_rng(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
    Y_pred = forwardsim(drug_regimen_pred,
      Hb_star+theta_rand_pred[1],
      alpha_diff1,
      alpha_delta1,
      alpha_diff2,
      alpha_delta2,
      logit_alpha+theta_rand_pred[2],
      h,
      beta*exp(theta_rand_pred[3]),
      T_E_star+theta_rand_pred[4],
      log_k,
      N_pred,
      T_nmblast, T_retic, T_RBC_max, T_transit_steady_state, 
      mean_delay,
      sigma_delay,
      K_weights,
      sigma);
  }
}


