require(gtools)

# This is an R version of the enzyme kinetics model - helps debugging
forwardsim = function(drug_regimen,              # dose given each day (mg/kg)
                      Hb_star,                   # Steady state hameoglobin for the individual (g/dL)
                      diff_alpha,                # parameters governing increase in normoblast production based on difference from steady state
                      delta_alpha,               # parameters governing increase in normoblast production based on previous day reduction in haemoglobin
                      MAX_EFFECT,                # maximum drug effect on G6PD depletion
                      h,                         # slope effect in dose response function
                      beta,                      # Dose corresponding to half maximum drug effect (mg)
                      log_k,                     # transit time parameter
                      nComp_sim,                 # Total number of compartments to forward simulate
                      T_nmblast,                 # normboblast lifespan - maturation time
                      T_retic,                   # reticulocyte lifespan - including marrow and circulating periods
                      T_RBC_max,                 # Maximum allowed value for RBC lifespan (arbitrary; days)
                      T_transit_steady_state,    # transit time at steady state
                      G6PD_initial,              # Initial amount of G6PD (unitless, scales with the delta quantity)
                      G6PD_decay_rate,      # Daily fixed reduction of G6PD in absence of drug
                      mu_death,
                      sigma_death
                      ){
  
  out_res = array(dim = c(2, nComp_sim))       # where we save the output of forwards simulation
  
  # set up variables we need for the simulation
  
  Hb_delta=0;
  lambda = 1000;                            # baseline production of normoblasts per time - arbitrary number
  
  Hb = array(dim=nComp_sim);               # The haemoglobin at each point in the simulation
  retic_percent = array(dim=nComp_sim)     # The retic count (%) at each point in the simulation
  drug_effect=0;
  
  # s to store red blood cell age distributions
  erythrocytes = rep(lambda, T_RBC_max);         # Distribution of erythrocytes in circulation
  normoblasts = rep(lambda, T_nmblast);          # Distribution of normoblasts
  reticulocytes = rep(lambda, T_retic);          # Distribution of reticulocytes
  erythrocytes_G6PD = array(dim = T_RBC_max)
  
  # ***** Calculate initial values at steady state *****
  Hb[1] = Hb_star;
  transit = compute_transit_time(Hb[1], Hb_star, T_transit_steady_state, log_k);
  rho = compute_rho(Hb[1], Hb_delta, Hb_star, diff_alpha, delta_alpha); # this should return 1
  if(rho>1) print("Problem as rho >1 at steady state!");
  
  erythrocytes_G6PD[1] = G6PD_initial; # starting quantity of G6PD enzyme
  for(i in 2:T_RBC_max){
    erythrocytes_G6PD[i] = erythrocytes_G6PD[i-1]*exp(-G6PD_decay_rate); # daily decay
    # killing probability based on G6PD enzyme quantity
    erythrocytes[i] = erythrocytes[i-1]*inv.logit((log(erythrocytes_G6PD[i-1])-mu_death)*sigma_death);
  }
  plot(erythrocytes_G6PD, type='l',lty=1,col='black',lwd=2)
  
  Total_Retics = CountRetics(transit, reticulocytes); # count the number of retics in circulation
  Total_Eryths = sum(erythrocytes);                   # count the number of mature RBCs in circulation
  C_0 = Total_Retics + Total_Eryths;                  # Starting number of circulating cells - this is the steady state number
  retic_percent[1] = 100 * Total_Retics / C_0;
  
  #******************* Forward simulation *****************#
  for(t in 2:nComp_sim){
    
    # Compute the multiplication factor of basal normoblast production
    if(t>2) Hb_delta = Hb[t-2]-Hb[t-1];
    rho = compute_rho(Hb[t-1], Hb_delta, Hb_star, diff_alpha, delta_alpha);
    
    # calculate the updated Hb dependent transit time for the reticulocytes
    transit = compute_transit_time(Hb[t-1], Hb_star, T_transit_steady_state, log_k);
    
    # We move the RBCs from one compartment to the next
    temp_normoblasts = normoblasts;   # make temp copy
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
    
    # Compute the G6PD depletion effect of yesterday's dose on circulating erythrocytes
    drug_effect = dose_response(drug_regimen[t-1], MAX_EFFECT, h, beta);
    
    # survival function of RBC lifespan
    temp_eryths = erythrocytes;             # copy values to temp 
    temp_eryths_G6PD = erythrocytes_G6PD;
    # values for first day as erythrocytes
    erythrocytes[1] = temp_retics[T_retic];
    erythrocytes_G6PD[1] = G6PD_initial;
    for (i in 2:T_RBC_max){                # move all the erythrocytes along by one, but now including G6PD enzyme quantity and induced killing rates
      # deplete G6PD enzyme quantity by daily amount
      erythrocytes_G6PD[i] = temp_eryths_G6PD[i-1]*exp(-G6PD_decay_rate-drug_effect); 
      # death probability determined by log enzyme quantity
      erythrocytes[i] = temp_eryths[i-1]*inv.logit((log(erythrocytes_G6PD[i-1])-mu_death)*sigma_death);
    }
    
    # Count the number of retics and erythrocytes in circulation
    Total_Retics = CountRetics(transit, reticulocytes);
    Total_Eryths = sum(erythrocytes);
    C_t = Total_Retics + Total_Eryths;
    retic_percent[t] = 100 * Total_Retics / C_t;
    Hb[t] = Hb_star * C_t / C_0;
    
  }
  # store results
  out_res[1,] = Hb;
  out_res[2,] = retic_percent;
  lines(erythrocytes_G6PD, type='l',lty=2,col='red')
  return(out_res);
}


CountRetics = function( transit,  reticulocytes){
  Total_Retics = 0;
  T_retic = length(reticulocytes);
  i = T_retic;
  transit_int = ceiling(transit);
  retic_rem =  ceiling(transit) - transit;
  while(i > transit+1) {
    Total_Retics = Total_Retics+reticulocytes[i];
    i = i-1;
  }
  Total_Retics = Total_Retics+retic_rem*reticulocytes[transit_int];
  
  return(Total_Retics);
}



compute_rho = function( Hb,  Hb_delta,  Hb_star,  diff_alpha,  delta_alpha)
{
  Hb_delta_temp=Hb_delta;
  Hb_diff_temp=Hb_star-Hb;
  rho=1; 
  K_diff = length(diff_alpha);
  K_delta = length(delta_alpha);
  
  if(Hb_delta_temp<0) Hb_delta_temp=0;
  for(i in 1:K_delta){
    rho =rho+ delta_alpha[i]*(Hb_delta_temp^i);
  }
  
  if(Hb_diff_temp<0) Hb_diff_temp=0;
  for(i in 1:K_diff){
    rho = rho + diff_alpha[i]*(Hb_diff_temp^i);
  }
  return (rho);
}

dose_response = function( dose,  MAX_EFFECT,  h,  beta){
  effect = MAX_EFFECT * (dose ^h ) / (dose ^ h + beta^h);
  return (effect);
}

compute_transit_time = function( Hb,  Hb_star,  T_transit_steady_state,  log_k)
{
  transit = T_transit_steady_state*exp(-(Hb_star-Hb)*exp(log_k));
  if(transit>T_transit_steady_state) transit=T_transit_steady_state;
  return (transit);
}

xs = seq(0, 1, length.out=100)
ys=dose_response(xs,
                 MAX_EFFECT = exp(-2), 
                 h = 4.5, 
                 beta = exp(-2.2))
plot(xs,ys,type='l')
lines(xs,dose_response(xs,
                      MAX_EFFECT = inv.logit(-2.5), 
                      h = 0.1, 
                      beta = exp(-1.6)))

dose_response(0.5,MAX_EFFECT = 0.01, h = 5, beta = 0.25)


plot((xs), ys,type='l')
nComp_sim = 100
out=forwardsim(drug_regimen = c(rep(0.5,nComp_sim)), 
               Hb_star = 15, 
               diff_alpha = c(.1), 
               delta_alpha = c(.5),
               MAX_EFFECT = exp(-2.5), 
               h = 5, 
               beta = 0.25,
               log_k = 1,
               nComp_sim = nComp_sim,
               T_nmblast = 5, T_retic = 5, T_RBC_max = 150, 
               T_transit_steady_state = 3.5,
               G6PD_initial = 1,
               G6PD_decay_rate = exp(-3),
               mu_death = -3, 
               sigma_death = 8)

par(mfrow=c(1,2))
plot(out[1,],type='l')
plot(out[2,],type='l')

require(gtools)
par(las=1, mfrow=c(1,2))
xs= 0:100
rate_decay = 0.05
plot(xs, exp(-rate_decay * xs), type='l', ylim = c(0,1),panel.first=grid(),
     xlab='Age (days)', ylab = 'Relative enzyme quantity')
plot(xs, 1-inv.logit((-rate_decay * xs + 2.5)*7.6), type='l',panel.first=grid(),
     ylab='Daily probability of death', xlab = 'Age (days)',ylim = c(0,1))
