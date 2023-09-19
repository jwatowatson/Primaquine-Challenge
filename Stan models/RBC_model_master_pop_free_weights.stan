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


// NOTE: Have replaced the parametric dose-weighting parameters with an
// independent `dose_weights` simplex.


functions {
  // This function computes the effective dose given a dosing schedule vector
  // drug_regimen: a vector of daily primaquine `equivalent' doses in mg
  // (in practice only dosed once a day, this can be changed to hours, weeks, ...)
  // t: current time point
  // NOTE: Have replaced the weighting parameters with the `weights` vector.
  vector compute_effective_dose(vector drug_regimen, int nComp_sim, vector weights) {
    vector[nComp_sim] effective_dose;

    for (t in 1:nComp_sim) {
      effective_dose[t] = 0;
    }

    for (t in 1:(nComp_sim-1)) {
      for (kk in 1:(nComp_sim-t)) {
        effective_dose[t + kk] += weights[kk] * drug_regimen[t];
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
    Total_retics += retic_rem*reticulocytes[transit_int];

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
  vector dose_weights,            // NOTE: weights used to calculate effective doses
  real mean_delay,
  real sigma_delay,
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
    // NOTE: Have replaced the weighting parameters with the `weights` vector.
    effective_dose = compute_effective_dose(drug_regimen, nComp_sim, dose_weights);

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
  real T_transit_steady_state;   // Time at which reticulocytes transit into circulation
  real sigma;                    // sd for cumulative standard normal distribution in RBC survival function

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
  real<lower=0> T_E_star_mean;
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
  // vector of zeros for the random effects specification
  int K_rand_effects = 8;
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
  real<lower=0> mean_delay;
  real<lower=1> sigma_delay;

  // NOTE: effective dose weights
  simplex[28] dose_weights;

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
  // Individual random effects
  vector[K_rand_effects] theta_rand[N_subject]; // individual random effects vector
  cholesky_factor_corr[K_rand_effects] L_Omega; // correlation matrix
  vector<lower=0>[K_rand_effects] sigmasq_u;    // variance of random effects

  // Repeated occasions random effects
  real theta_ic[N_repeat];       // repeated occasions random effect
  real<lower=0> sigmasq_u_ic;    // variance of random effects inter-occasion
}

transformed parameters {
  matrix[2,N_sim_tot] Y_hat;

  for(j in 1:N_experiment){
    vector[K_rand_effects] theta_ic_j = theta_rand[id[j]];
    if(repeat_oc[j]>0){
      // then we add inter-occasion variability to the following parameters:
      // the baseline haemoglobin (steady state?)
      theta_ic_j[1] += theta_ic[repeat_oc[j]];
    }
    // compute output from forward_sim
    Y_hat[,ind_start_regimen[j]:ind_end_regimen[j]] = forwardsim(
      drug_regimen[ind_start_regimen[j]:ind_end_regimen[j]], // drug doses
      Hb_star + theta_ic_j[1],                               // steady state haemoglobin
      alpha_diff1*exp(theta_ic_j[2]),                        // parameters on bone marrow response (polynomial)
      alpha_delta1*exp(theta_ic_j[3]),
      alpha_diff2*exp(theta_ic_j[4]),
      alpha_delta2*exp(theta_ic_j[5]),
      logit_alpha+theta_ic_j[6],                             // max effect on lifespan
      h,                                                     // slope of effect on lifespan
      beta*exp(theta_ic_j[7]),                               // dose giving half max effect on lifespan
      T_E_star + theta_ic_j[8],                              // steady state lifespan
      log_k,                                                 // retic release parameter
      N_sim[j],
      T_nmblast,
      T_retic,
      T_RBC_max,
      T_transit_steady_state,
      dose_weights,
      mean_delay,
      sigma_delay,
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
  sigmasq_u[1] ~ normal(1,1) T[0,];    // Hb_star
  sigmasq_u[2] ~ exponential(10);      // alpha_diff1
  sigmasq_u[3] ~ exponential(10);      // alpha_delta1
  sigmasq_u[4] ~ exponential(10);      // alpha_diff2
  sigmasq_u[5] ~ exponential(10);      // alpha_delta2
  sigmasq_u[6] ~ normal(0.5,.5) T[0,]; // logit_alpha
  sigmasq_u[7] ~ exponential(10);      // beta
  sigmasq_u[8] ~ exponential(1);       // T_E_star

  sigmasq_u_ic ~ normal(.5, .5);

  L_Omega ~ lkj_corr_cholesky(2);
  for(i in 1:N_subject) theta_rand[i] ~ multi_normal_cholesky(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
  theta_ic ~ normal(0, sigmasq_u_ic);

  // steady state parameters
  T_E_star ~ normal(T_E_star_mean,T_E_star_sigma);
  Hb_star ~ normal(Hb_star_mean,Hb_star_sigma);

  // delay in dose effect
  mean_delay ~ exponential(.1);
  sigma_delay ~ exponential(.1);

  // NOTE: effective dose weights
  dose_weights ~ dirichlet(rep_vector(1, 28));

  // parameters governing the dose-response curve
  h ~ exponential(1);
  logit_alpha ~ normal(alpha_mean, alpha_sigma);
  beta ~ normal(beta_mean, beta_sigma) T[0,];

  // parameters governing the production of new cells in bone marrow
  alpha_delta1 ~ exponential(5);
  alpha_diff1 ~  exponential(5);
  alpha_delta2 ~ exponential(5);
  alpha_diff2 ~  exponential(5);

  // retic transit function
  log_k ~ normal(log_k_mean,log_k_sigma);

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

    theta_rand_pred = multi_normal_cholesky_rng(my_zeros, diag_pre_multiply(sigmasq_u, L_Omega));
    Y_pred = forwardsim(drug_regimen_pred,
    Hb_star+theta_rand_pred[1],
    alpha_diff1*exp(theta_rand_pred[2]),
    alpha_delta1*exp(theta_rand_pred[3]),
    alpha_diff2*exp(theta_rand_pred[4]),
    alpha_delta2*exp(theta_rand_pred[5]),
    logit_alpha+theta_rand_pred[6],
    h,
    beta*exp(theta_rand_pred[7]),
    T_E_star+theta_rand_pred[8],
    log_k,
    N_pred,
    T_nmblast,
    T_retic,
    T_RBC_max,
    T_transit_steady_state,
    dose_weights,
    mean_delay,
    sigma_delay,
    sigma);
  }
}
