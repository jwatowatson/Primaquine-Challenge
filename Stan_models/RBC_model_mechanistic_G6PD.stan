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
    effect = MAX_EFFECT * pow(dose,h)/(pow(dose,h) + pow(beta,h)); // EMAX function
    return effect;
  }
  
  // Compute the fold change from baseline in the number of new normoblasts produced given the Hb level
  // Hb: Hb level at time t
  // Hb_star: steady-state Hb
  // rho_max: maximum fold increase in production relative to steady state
  real compute_rho(real Hb, real Hb_delta, real Hb_star, real diff_alpha, real delta_alpha)
  {
    real rho;
    rho = exp(diff_alpha*(Hb-Hb_star))*exp(delta_alpha*Hb_delta);
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
      i += -1;
    }
    Total_Retics += retic_rem*reticulocytes[transit_int];
    
    return Total_Retics;
  }
  
  // ************ deterministic simulation using the RBC dynamics model parameters **************
  matrix forwardsim(
    vector drug_regimen,            // dose given each day (mg/kg)
    real Hb_star,                   // Steady state hameoglobin for the individual (g/dL)
    real diff_alpha,              // parameters governing increase in normoblast production based on difference from steady state
    real delta_alpha,             // parameters governing increase in normoblast production based on previous day reduction in haemoglobin
    real MAX_EFFECT,                // maximum additional rate of decay in G6PD from drug effect 
    real h,                         // slope effect in dose response function
    real beta,                      // Dose corresponding to half maximum drug effect (mg)
    real log_k,                     // transit time parameter
    int  nComp_sim,                 // Total number of compartments to forward simulate
    int  T_nmblast,                 // normboblast lifespan - maturation time
    int  T_retic,                   // reticulocyte lifespan - including marrow and circulating periods
    int  T_RBC_max,                 // Maximum allowed value for RBC lifespan (arbitrary; days)
    real T_transit_steady_state,    // transit time at steady state
    real G6PD_initial,              // Initial amount of G6PD (unitless, scales with the delta quantity)
    real G6PD_decay_rate,            // Daily fixed proportion of G6PD depleted in absence of drug
    real mu_death,
    real sigma_death,
    int K_out
    ){
      
      matrix[K_out, nComp_sim] out_res;            // for saving the output of the forwards simulation
      
      // set up variables we need for the simulation
      real transit;
      row_vector[nComp_sim]  rho;              // increase in production in bone marrow
      real Hb_delta=0;
      real lambda = 1000;                      // baseline production of normoblasts per time - arbitrary number
      real Total_Eryths;                       // Total circulating erythrocytes
      real Total_Retics;                       // Total circulating reticulocytes
      
      row_vector[nComp_sim] Hb;                // The haemoglobin at each point in the simulation
      row_vector[nComp_sim] retic_percent;     // The retic count (%) at each point in the simulation
      row_vector[nComp_sim] drug_effect;
      real C_t;                                // Counter of cells
      real C_0;                                // Counter of cells - at baseline steady state
      
      // vectors to store red blood cell age distributions
      vector[T_RBC_max] erythrocytes = rep_vector(lambda, T_RBC_max);         // Distribution of erythrocytes in circulation
      vector[T_RBC_max] temp_eryths; // need the temp vectors to store temp info otherwise it gets overwritten
      vector[T_nmblast] normoblasts = rep_vector(lambda, T_nmblast);          // Distribution of normoblasts
      vector[T_nmblast] temp_normoblasts;
      vector[T_retic] reticulocytes = rep_vector(lambda, T_retic);            // Distribution of reticulocytes
      vector[T_retic] temp_retics;
      
      // vectors to store G6PD total amount per erythrocyte age
      vector[T_RBC_max] erythrocytes_G6PD;    // G6PD in erythrocytes in circulation
      vector[T_RBC_max] temp_eryths_G6PD;
      
      
      // ***** Calculate initial values at steady state *****
      Hb[1] = Hb_star;
      transit = compute_transit_time(Hb[1], Hb_star, T_transit_steady_state, log_k);
      rho[1] = compute_rho(Hb[1], Hb_delta, Hb_star, diff_alpha, delta_alpha); // this should return 1 - this value is not used
      if(rho[1]>1) print("*******Problem as rho >1 at steady state!*********");
      drug_effect[1]=0.0;
      
      erythrocytes_G6PD[1] = G6PD_initial; // starting quantity of G6PD enzyme
      for(i in 2:T_RBC_max){
        erythrocytes_G6PD[i] = erythrocytes_G6PD[i-1]*exp(-G6PD_decay_rate); // daily decay
        // killing probability based on G6PD enzyme quantity
        erythrocytes[i] = erythrocytes[i-1]*inv_logit((log(erythrocytes_G6PD[i-1])-mu_death)*sigma_death);
      }
      
      Total_Retics = CountRetics(transit, reticulocytes); // count the number of retics in circulation
      Total_Eryths = sum(erythrocytes);                   // count the number of mature RBCs in circulation
      C_0 = Total_Retics + Total_Eryths;                  // Starting number of circulating cells - this is the steady state number
      retic_percent[1] = 100 * Total_Retics / C_0;
      
      //******************* Forward simulation *****************//
      for(t in 2:nComp_sim){
        
        // Compute the multiplication factor of basal normoblast production
        if(t>2) Hb_delta = Hb[t-2]-Hb[t-1];
        rho[t] = compute_rho(Hb[t-1], Hb_delta, Hb_star, diff_alpha, delta_alpha);
        
        // calculate the updated Hb dependent transit time for the reticulocytes
        transit = compute_transit_time(Hb[t-1], Hb_star, T_transit_steady_state, log_k);
        
        // We move the red cells from one compartment to the next
        temp_normoblasts = normoblasts;   // make temp copy
        normoblasts[1] = rho[t]*lambda;   // the number of new normoblasts made at time t
        for(i in 2:T_nmblast)             // move all the normoblasts along by one
        {
          normoblasts[i] = temp_normoblasts[i-1];
        }
        
        temp_retics = reticulocytes;
        reticulocytes[1] = temp_normoblasts[T_nmblast];
        for(i in 2:T_retic)              // move all the retics along by one
        {
          reticulocytes[i] = temp_retics[i-1];
        }
        
        // Compute the additional G6PD decay effect from yesterday's primaquine dose - EMAX function
        drug_effect[t] = dose_response(drug_regimen[t-1], MAX_EFFECT, h, beta);
        
        // survival function of RBC lifespan
        temp_eryths = erythrocytes; // copy values to temp vector
        temp_eryths_G6PD = erythrocytes_G6PD;
        // values for first day as erythrocytes
        erythrocytes[1] = temp_retics[T_retic];
        erythrocytes_G6PD[1] = G6PD_initial;
        for (i in 2:T_RBC_max){ // move all the erythrocytes along by one, but now including G6PD enzyme quantity and induced killing rates
          // deplete G6PD enzyme quantity by daily amount
          erythrocytes_G6PD[i] = temp_eryths_G6PD[i-1]*exp(-G6PD_decay_rate-drug_effect[t]); 
          // death probability determined by log enzyme quantity
          erythrocytes[i] = temp_eryths[i-1]*inv_logit((log(erythrocytes_G6PD[i-1])-mu_death)*sigma_death);
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
      out_res[3] = rho;
      out_res[4] = drug_effect;
      
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
  real log_beta_mean;
  real<lower=0> Hb_star_mean;
  real log_k_mean;
  real log_MAX_EFFECT_prior_mean;
  
  //variation
  real<lower=0> log_slope_effect_var;
  real<lower=0> beta_sigma;
  real<lower=0> Hb_star_sigma;
  real<lower=0> log_k_sigma;
  real<lower=0> log_MAX_EFFECT_prior_sigma;
  
  //residual
  real<lower=0> sigma_Hb_mean;
  
  // single drug regimen to predict
  int<lower=1> N_pred;
  vector[N_pred] drug_regimen_pred;
  
}

transformed data{
  // vector of zeros for the random effects specification
  int K_rand_effects = 6;
  real G6PD_initial=1.0;
  real mu_death=-3.0;             // governing death process
  int K_out=4;
  vector[K_rand_effects] my_zeros;
  real sigma_retic=1;
  for(i in 1:K_rand_effects) my_zeros[i] = 0;
}

parameters {
  
  // error terms
  real<lower=0> sigma_CBC;
  real<lower=0> sigma_haemocue;
  // real<lower=0> sigma_retic;
  real CBC_correction; // systematic difference between CBC haemoglobin and haemocue haemoglobin
  
  // parameters governing the dose-response curve
  real log_MAX_EFFECT;   // maximal effect on log scale (max log % depletion)
  real log_beta;         // half maximal effect on the log mg/kg scale
  real<lower=1> h;       // hill coefficient for EMAX function
  
  // parameter governing G6PD depletion (starting amount versus daily reduction are not identifiable so we fix daily depletion at 1)
  real log_G6PD_decay_rate;  // population mean daily log decay rate
  real<lower=0> sigma_death; // governing death process
  
  // retic transit function
  real log_k;
  
  real<lower=0> Hb_star; // steady state haemoglobin
  
  // parameters governing the production of new cells in bone marrow
  real diff_alpha;
  real delta_alpha;
  
  // random effects
  // Individual random effects
  vector[K_rand_effects] theta_rand[N_subject]; // individual random effects vector
  cholesky_factor_corr[K_rand_effects] L_Omega; // correlation matrix
  vector<lower=0>[K_rand_effects] sigmasq_u;    // variance of random effects
  
}

transformed parameters {
  matrix[K_out,N_sim_tot] Y_hat;
  for(j in 1:N_experiment){
    // compute output from forward_sim
    Y_hat[,ind_start_regimen[j]:ind_end_regimen[j]] = 
    forwardsim(
      drug_regimen[ind_start_regimen[j]:ind_end_regimen[j]], // drug doses
      Hb_star + theta_rand[id[j]][1],                        // steady state haemoglobin
      diff_alpha+theta_rand[id[j]][4],                  // bone marrow response as function of difference from steady state Hb
      delta_alpha+theta_rand[id[j]][5],                 // bone marrow response as function of daily change in Hb
      exp(log_MAX_EFFECT+theta_rand[id[j]][2]),              // max effect on G6PD decay rate
      h,                                                     // slope of effect on lifespan
      exp(log_beta+theta_rand[id[j]][3]),                    // dose giving half max effect on lifespan
      log_k,                                                 // retic release parameter
      N_sim[j],               // Number of days to forward simulate
      T_nmblast,              // days as normoblasts: FIXED=5
      T_retic,                // days as reticulocytes: FIXED=5
      T_RBC_max,              // Max length of erythrocyte vector: FIXED=150
      T_transit_steady_state, // Time of transit into circulation: FIXED=3.5
      G6PD_initial,           // FIXED=1            
      exp(log_G6PD_decay_rate+theta_rand[id[j]][6]),// daily rate of decay of G6PD enzyme
      mu_death,               // FIXED=-3 log enzyme quantity resulting in daily death rate of 50%
      sigma_death,            // slope coefficient around the 50% mark
      K_out);
  }
}

model{
  // Prior
  // error terms
  sigma_CBC ~ exponential(1);
  sigma_haemocue ~ exponential(1);
  // sigma_retic ~ normal(0.5, .25) T[0,]; // VERY STRONG PRIOR-DO WE NEED THIS?
  CBC_correction ~ normal(0,1);
  
  // random effects
  sigmasq_u[1] ~ normal(1,1) T[0,];    // Hb_star
  sigmasq_u[2] ~ exponential(20);      // MAX_EFFECT
  sigmasq_u[3] ~ exponential(20);      // beta
  sigmasq_u[4] ~ exponential(1);       // alpha_diff
  sigmasq_u[5] ~ exponential(1);       // alpha_delta
  sigmasq_u[6] ~ exponential(10);      // individual G6PD decay rate in absence of drug
  
  L_Omega ~ lkj_corr_cholesky(2);
  for(i in 1:N_subject) theta_rand[i] ~ multi_normal_cholesky(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));

  // steady state parameters
  Hb_star ~ normal(Hb_star_mean,Hb_star_sigma);
  
  // parameters governing the dose-response curve
  h ~ normal(4,1) T[1,];
  log_MAX_EFFECT ~ normal(log_MAX_EFFECT_prior_mean, log_MAX_EFFECT_prior_sigma);
  log_beta ~ normal(log_beta_mean, beta_sigma);
  
  // parameters governing the production of new cells in bone marrow
  delta_alpha ~ normal(0,1);
  diff_alpha ~ normal(0,1);
  
  // retic transit function
  log_k ~ normal(log_k_mean,log_k_sigma);
  
  // G6PD depletion process - daily decay rate process
  log_G6PD_decay_rate ~ normal(-4.6,1); // prior is 1% per day 
  sigma_death ~ normal(10, 5);
  
  // Likelihood
  for(j in 1:N_experiment){
    // calculate likelihood
    Hb_CBC[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]] ~ normal(Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_CBC_Hb[ind_start_Hb_CBC[j]:ind_end_Hb_CBC[j]]]+CBC_correction, sigma_CBC);
    Hb_Haemocue[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]] ~ normal(Y_hat[1][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_hemocue_Hb[ind_start_Hb_hemocue[j]:ind_end_Hb_hemocue[j]]], sigma_haemocue);
    // POSSIBLY NOT BEST LIKELIHOOD FOR RETICS
    Retic_data[ind_start_retic[j]:ind_end_retic[j]] ~ normal(Y_hat[2][ind_start_regimen[j]:ind_end_regimen[j]][t_sim_retic[ind_start_retic[j]:ind_end_retic[j]]], sigma_retic);
  }
  
}


generated quantities {
  
  matrix[K_out,N_pred] Y_pred;
  {
    vector[K_rand_effects] theta_rand_pred; // individual random effects vector
    
    theta_rand_pred = multi_normal_cholesky_rng(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
    Y_pred = forwardsim(
      drug_regimen_pred,                     // drug doses
      Hb_star + theta_rand_pred[1],          // steady state haemoglobin
      diff_alpha+theta_rand_pred[4],         // parameters on bone marrow response (polynomial)
      delta_alpha+theta_rand_pred[5],
      exp(log_MAX_EFFECT+theta_rand_pred[2]),// max effect on lifespan
      h,                                     // slope of effect on lifespan
      exp(log_beta+theta_rand_pred[3]),      // dose giving half max effect on lifespan
      log_k,                                 // retic release parameter
      N_pred,                                //FIXED
      T_nmblast,                             //FIXED
      T_retic,                               //FIXED
      T_RBC_max,                             //FIXED
      T_transit_steady_state,                //FIXED 
      G6PD_initial,                          //FIXED           
      exp(log_G6PD_decay_rate+theta_rand_pred[6]),
      mu_death,                              //FIXED
      sigma_death,
      K_out);
  }
}
