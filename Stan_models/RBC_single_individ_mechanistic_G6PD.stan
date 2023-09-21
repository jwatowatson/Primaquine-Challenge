//
//###############################################################################
//# Authors: James Watson and Rob Moss                                                 #
//# Date: Sept 2023                                                                    #
//# Name: RBC mechanistic                                                              #
//# Description: Mechanistic version of model explicitly modelling G6PD depletion      #
//###############################################################################
//# This code includes two parts including the red blood cell dynamics functions
//# and the hierarchial model fitting


functions {
  
  //function for changing real to integer used for calculating transit time
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
  
  
  // This function computes the reduction in RBC lifespan
  // MAX_EFFECT: maximum % decrease in RBC lifespan - logit scale
  // h: slope of the curve in log scale
  // beta: effective dose value when maximum drug effect is halved.
  real dose_response(real dose, real MAX_EFFECT, real h, real beta){
    real effect;
    effect = MAX_EFFECT * pow(dose,h)/ (pow(dose,h) + pow(beta,h)); // EMAX function
    return effect;
  }
  
  // Compute the fold change from baseline in the number of new normoblasts produced given the Hb level
  // Hb: Hb level at time t
  // Hb_star: steady-state Hb
  // rho_max: maximum fold increase in production relative to steady state
  // Hb_rho_half: mid-point Hb when maximum production is halved
  real compute_rho(real Hb, real Hb_delta, real Hb_star, vector diff_alpha, vector delta_alpha)
  {
    real Hb_delta_temp=Hb_delta;
    real Hb_diff_temp=Hb_star-Hb;
    real rho=1; // rho = 1 is baseline production
    int K_diff = num_elements(diff_alpha);
    int K_delta = num_elements(delta_alpha);
    
    if(Hb_delta_temp<0) Hb_delta_temp=0;
    for(i in 1:K_delta){
      rho += delta_alpha[i]*pow(Hb_delta_temp,i);
    }
    
    if(Hb_diff_temp<0) Hb_diff_temp=0;
    for(i in 1:K_diff){
      rho += diff_alpha[i]*pow(Hb_diff_temp, i);
    }
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
  
  
  // Count the number of retics in circulation
  // transit: transit time
  // reticulocytes: retics vector at time t
  real CountRetics(real transit, vector reticulocytes){
    real Total_Retics = 0;
    int T_retic = num_elements(reticulocytes);
    int i = T_retic;
    real transit_idx = ceil(transit);
    int transit_int = real2int(transit_idx, 0, 6);
    real retic_rem =  ceil(transit) - transit;
    while(i > transit+1) {
      Total_Retics += reticulocytes[i];
      i = i-1;
    }
    Total_Retics += retic_rem*reticulocytes[transit_int];
    
    return Total_Retics;
  }
  
  // ************ deterministic simulation using the RBC dynamics model parameters **************
  matrix forwardsim(
    vector drug_regimen,            // dose given each day (mg/kg)
    real Hb_star,                   // Steady state hameoglobin for the individual (g/dL)
    vector diff_alpha,              // parameters governing increase in normoblast production based on difference from steady state
    vector delta_alpha,             // parameters governing increase in normoblast production based on previous day reduction in haemoglobin
    real MAX_EFFECT,                // maximum drug effect on G6PD depletion
    real h,                         // slope effect in dose response function
    real beta,                      // Dose corresponding to half maximum drug effect (mg)
    real log_k,                     // transit time parameter
    int  nComp_sim,                 // Total number of compartments to forward simulate
    int  T_nmblast,                 // normboblast lifespan - maturation time
    int  T_retic,                   // reticulocyte lifespan - including marrow and circulating periods
    int  T_RBC_max,                 // Maximum allowed value for RBC lifespan (arbitrary; days)
    real T_transit_steady_state,    // transit time at steady state
    real G6PD_initial,              // Initial amount of G6PD (unitless, scales with the delta quantity)
    real G6PD_delta_day,      // Daily fixed proportion of G6PD depleted in absence of drug
    real G6PD_threshold
    ){
      
      matrix[2, nComp_sim] out_res;           // where we save the output of forwards simulation
      
      // set up variables we need for the simulation
      real transit;
      real rho;
      real Hb_delta=0;
      real lambda = 100;                      // baseline production of normoblasts per time - arbitrary number
      real Total_Eryths;                      // Total circulating erythrocytes
      real Total_Retics;                      // Total circulating reticulocytes
      
      row_vector[nComp_sim] Hb;               // The haemoglobin at each point in the simulation
      row_vector[nComp_sim] retic_percent;    // The retic count (%) at each point in the simulation
      real drug_effect=0;
      real C_t;                               // Counter of cells
      real C_0;                               // Counter of cells - at baseline steady state
      
      // vectors to store red blood cell age distributions
      vector[T_RBC_max] erythrocytes = rep_vector(lambda, T_RBC_max);         // Distribution of erythrocytes in circulation
      vector[T_RBC_max] temp_eryths; // need the temp vectors to store temp info otherwise it gets overwritten
      vector[T_nmblast] normoblasts = rep_vector(lambda, T_nmblast);          // Distribution of normoblasts
      vector[T_nmblast] temp_normoblasts;
      vector[T_retic] reticulocytes = rep_vector(lambda, T_retic);          // Distribution of reticulocytes
      vector[T_retic] temp_retics;
      
      // vectors to store G6PD total amount per erythrocyte age
      vector[T_RBC_max] erythrocytes_G6PD;    // G6PD in erythrocytes in circulation
      vector[T_RBC_max] temp_eryths_G6PD;
      
      
      // ***** Calculate initial values at steady state *****
      Hb[1] = Hb_star;
      transit = compute_transit_time(Hb[1], Hb_star, T_transit_steady_state, log_k);
      rho = compute_rho(Hb[1], Hb_delta, Hb_star, diff_alpha, delta_alpha); // this should return 1
      if(rho>1) print("Problem as rho >1 at steady state!");
      
      erythrocytes_G6PD[1] = G6PD_initial; // starting quantity of G6PD enzyme
      for(i in 2:T_RBC_max){
        erythrocytes_G6PD[i] = erythrocytes_G6PD[i-1]*(1-G6PD_delta_day); // daily depletion
        if(erythrocytes_G6PD[i-1] < G6PD_threshold){
          erythrocytes[i] = 0.0;
        } else {
          erythrocytes[i] = erythrocytes[i-1];
        }
      }
      
      Total_Retics = CountRetics(transit, reticulocytes); // count the number of retics in circulation
      Total_Eryths = sum(erythrocytes);                   // count the number of mature RBCs in circulation
      C_0 = Total_Retics + Total_Eryths;                  // Starting number of circulating cells - this is the steady state number
      retic_percent[1] = 100 * Total_Retics / C_0;
      
      //******************* Forward simulation *****************//
      for(t in 2:nComp_sim){
        
        // Compute the multiplication factor of basal normoblast production
        if(t>2) Hb_delta = Hb[t-2]-Hb[t-1];
        rho = compute_rho(Hb[t-1], Hb_delta, Hb_star, diff_alpha, delta_alpha);
        
        // calculate the updated Hb dependent transit time for the reticulocytes
        transit = compute_transit_time(Hb[t-1], Hb_star, T_transit_steady_state, log_k);
        
        // We move the RBCs from one compartment to the next
        temp_normoblasts = normoblasts;   // make temp copy
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
        
        // Compute the G6PD depletion effect of yesterday's dose on circulating erythrocytes
        drug_effect = dose_response(drug_regimen[t-1], MAX_EFFECT, h, beta);
        
        // survival function of RBC lifespan
        temp_eryths = erythrocytes;             // copy values to temp vector
        temp_eryths_G6PD = erythrocytes_G6PD;
        // values for first day as erythrocytes
        erythrocytes[1] = temp_retics[T_retic];
        erythrocytes_G6PD[1] = G6PD_initial;
        for (i in 2:T_RBC_max){
          erythrocytes_G6PD[i] = temp_eryths_G6PD[i-1]*(1-G6PD_delta_day) - drug_effect; // deplete by daily amount plus drug effect
          if(erythrocytes_G6PD[i-1] < G6PD_threshold){
            erythrocytes[i] = 0.0;
          } else {
            erythrocytes[i] = temp_eryths[i-1];
          }
        }
        
        // Count the number of retics and erythrocytes in circulation
        Total_Retics = CountRetics(transit, reticulocytes);
        Total_Eryths = sum(erythrocytes);
        C_t = Total_Retics + Total_Eryths;
        retic_percent[t] = 100 * Total_Retics / C_t;
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
  real T_transit_steady_state;   // Time at which reticulocytes transit into circulation
  
  // Defining data to be fitted
  // by "experiment" we mean an individual trial of a single or ascending dose of primaquine
  // some volunteers have data both on single dose and on ascending dose - two seperate experiments
  int<lower=1> N_experiment;     // number of experiments (each subject can have 1 or more)
  int<lower=1,upper=N_experiment> N_subject; // number of unique subjects
  int<lower=0> N_repeat;        // number of repeated occasions
  int<lower=0> N_sim_tot;       // total number of days to simulate across all patients
  int<lower=0> N_CBC;           // total number of Hb - CBC observations
  int<lower=0> N_hemo;          // total number of Hb - Haemocue observations
  int<lower=0> N_retic;         // total number of retic (average of CBC and manual readings)
  
  // observed data and indexing
  // drug regimen to use and indexing
  int<lower=1> N_sim[N_experiment];   // length of simulation for each experiment
  int<lower=1> ind_start_regimen[N_experiment];  // starting indices for each experiment
  int<lower=1> ind_end_regimen[N_experiment];    // ending indices for each experiment
  vector[N_sim_tot] drug_regimen;                          // doses array (all individuals)
  int<lower=1,upper=N_subject> id[N_experiment];           // vector for subject ID (used for random effect terms)
  int<lower=0,upper=N_repeat> repeat_oc[N_experiment];     // vector for repeated occasions
  
  // - Hb CBC
  int<lower=1> ind_start_Hb_CBC[N_experiment];
  int<lower=1> ind_end_Hb_CBC[N_experiment];
  int<lower=1> t_sim_CBC_Hb[N_CBC];   // days where CBC haemoglobin data were available
  real<lower=0> Hb_CBC[N_CBC];        // CBC haemoglobin data
  
  // - Hb manual read
  int<lower=1> ind_start_Hb_hemocue[N_experiment];
  int<lower=1> ind_end_Hb_hemocue[N_experiment];
  int<lower=1> t_sim_hemocue_Hb[N_hemo];   // days where HaemoCue haemoglobin data were available
  real<lower=0> Hb_Haemocue[N_hemo];       // HaemoCue haemoglobin data
  
  // - Reticulocyte data
  int<lower=1> ind_start_retic[N_experiment];
  int<lower=1> ind_end_retic[N_experiment];
  int<lower=1> t_sim_retic[N_retic];       // days where reticulocyte data were available
  real<lower=0> Retic_data[N_retic];       // reticulocyte data
  
  // Prior parameters
  //mean
  real beta_mean;
  real<lower=0> Hb_star_mean;
  real log_k_mean;
  real logit_MAX_EFFECT_prior_mean;
  
  //variation
  real<lower=0> log_slope_effect_var;
  real<lower=0> beta_sigma;
  real<lower=0> Hb_star_sigma;
  real<lower=0> log_k_sigma;
  real<lower=0> logit_MAX_EFFECT_prior_sigma;
  
  //residual
  real<lower=0> sigma_Hb_mean;
  
  // single drug regimen to predict
  int<lower=1> N_pred;
  vector[N_pred] drug_regimen_pred;
  
}

transformed data{
  // vector of zeros for the random effects specification
  int K_rand_effects = 8;
  real G6PD_initial=1.0;
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
  real logit_MAX_EFFECT;   // maximal effect on logit scale (max logit % depletion)
  real<lower=0> beta;               // half maximal effect on the mg/kg scale
  real<lower=0> h;                  // hill coefficient for EMAX function
  
  // parameter governing G6PD depletion (starting amount versus daily reduction are not identifiable so we fix daily depletion at 1)
  real logit_G6PD_delta_day; // we fix the daily depletion amount - it is completely co-linear at steady state with the starting quantity
  real logit_G6PD_threshold; // governing death process
  
  // retic transit function
  real log_k;
  
  real<lower=0> Hb_star; // steady state haemoglobin
  
  // parameters governing the production of new cells in bone marrow
  vector<lower=0>[1] diff_alpha;
  vector<lower=0>[1] delta_alpha;
  
}

transformed parameters {
  matrix[2,N_sim_tot] Y_hat;
  
  for(j in 1:N_experiment){
    
    // compute output from forward_sim
    Y_hat[,ind_start_regimen[j]:ind_end_regimen[j]] = 
    forwardsim(
      drug_regimen[ind_start_regimen[j]:ind_end_regimen[j]], // drug doses
      Hb_star,                          // steady state haemoglobin
      diff_alpha,                       // parameters on bone marrow response (polynomial)
      delta_alpha,
      inv_logit(logit_MAX_EFFECT),      // max effect on lifespan
      h,                                // slope of effect on lifespan
      beta,                             // dose giving half max effect on lifespan
      log_k,                            // retic release parameter
      N_sim[j],
      T_nmblast,
      T_retic,
      T_RBC_max,
      T_transit_steady_state,
      G6PD_initial,              
      inv_logit(logit_G6PD_delta_day),
      inv_logit(logit_G6PD_threshold)
      );
  }
}

model{
  // Prior
  // error terms
  sigma_CBC ~ exponential(1);
  sigma_haemocue ~ exponential(1);
  sigma_retic ~ normal(0.5, .5) T[0,]; // VERY STRONG PRIOR-DO WE NEED THIS?
  CBC_correction ~ normal(0,1);

  // steady state parameters
  Hb_star ~ normal(Hb_star_mean,Hb_star_sigma);
  
  // parameters governing the dose-response curve
  h ~ normal(2, 0.5);
  logit_MAX_EFFECT ~ normal(logit_MAX_EFFECT_prior_mean, logit_MAX_EFFECT_prior_sigma);
  beta ~ normal(beta_mean, beta_sigma) T[0,];
  
  // parameters governing the production of new cells in bone marrow
  delta_alpha ~ exponential(10);
  diff_alpha ~  exponential(10);
  
  // retic transit function
  log_k ~ normal(log_k_mean,log_k_sigma);
  
  // G6PD depletion process - logit percent per day
  logit_G6PD_delta_day ~ normal(-3,0.25); // 
  logit_G6PD_threshold ~ normal(-5, 0.25);
  
  // Likelihood
  for(j in 1:N_experiment){
    // calculate likelihood
    Hb_CBC[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]] ~ student_t(10,Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_CBC_Hb[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]]]+CBC_correction, sigma_CBC);
    Hb_Haemocue[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]] ~ student_t(10,Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_hemocue_Hb[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]]], sigma_haemocue);
    log(Retic_data[ind_start_retic[j]:ind_end_retic[j]]) ~ normal(log(Y_hat[2][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_retic[ind_start_retic[j]:ind_end_retic[j]]]), sigma_retic);
  }
  
}


generated quantities {
  
  matrix[2,N_pred] Y_pred;
  {
    vector[K_rand_effects] theta_rand_pred; // individual random effects vector
    
    Y_pred = forwardsim(
      drug_regimen_pred, // drug doses
      Hb_star ,          // steady state haemoglobin
      diff_alpha,                            // parameters on bone marrow response (polynomial)
      delta_alpha,
      inv_logit(logit_MAX_EFFECT),        // max effect on lifespan
      h,                                     // slope of effect on lifespan
      beta,          // dose giving half max effect on lifespan
      log_k,                                 // retic release parameter
      N_pred,
      T_nmblast,
      T_retic,
      T_RBC_max,
      T_transit_steady_state,
      G6PD_initial,              
      inv_logit(logit_G6PD_delta_day),
      inv_logit(logit_G6PD_threshold)
      );
  }
}
