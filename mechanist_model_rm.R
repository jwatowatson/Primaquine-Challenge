#!/usr/bin/env -S Rscript --vanilla
#
# R implementation of `Stan_models/RBC_model_mechanistic_G6PD_rm.stan`.
#


main <- function() {
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(ggplot2)

  plot_file <- "mechanist_model_rm.png"

  job_id <- 6  # Single individual
  job_data <- load_job_data(job_id)
  drug_regimen <- c(extract_drug_regimen(job_data), rep(0, 30))
  clin_data <- extract_clinical_data(job_data)

  clin_all <- extract_all_data()

  # NOTE:
  # rho <- exp(diff_alpha * (Hb - Hb_star)) * exp(delta_alpha * Hb_delta)

  results <- forwardsim(
    drug_regimen = drug_regimen,
    Hb_star = 15,
    # diff_alpha = 0.1,
    # ## delta_alpha = 0.5,
    # ## delta_alpha = 3.5,
    # delta_alpha = 3.75,
    diff_alpha = 0.02,
    delta_alpha = 3.3,
    # MAX_EFFECT = exp(-2.5),
    MAX_EFFECT = exp(-2.75),
    ## h = 4,
    h = 3,
    ## beta = exp(-2),
    beta = exp(-1.8),
    log_k = -1,
    nComp_sim = length(drug_regimen),
    T_nmblast = 5,
    T_retic = 5,
    T_RBC_max = 120,
    T_transit_steady_state = 3.5,
    G6PD_initial = 1,
    G6PD_decay_rate = 0.01,
    mu_death = log(0.31),
    sigma_death = 70
  )

  p <- ggplot(
    results |> pivot_longer(! time),
    aes(time, value)
  ) +
    geom_line(
      aes(colour = factor(experiment)),
      data = clin_all,
      alpha = 0.4,
      show.legend = FALSE) +
    geom_line() +
    # geom_point(, data = clin_data$all) +
    expand_limits(y = 0) +
    facet_wrap(~ name, scales = "free_y", ncol = 1)

  cat("Writing", plot_file, "...\n")
  ggsave(plot_file, p, width = 15, height = 25, unit = "cm")

  plot_file <- "mechanist_model_rm_pred.png"
  drug_regimen <- c(job_data$drug_regimen_pred, rep(0, 30))

  results <- forwardsim(
    drug_regimen = drug_regimen,
    Hb_star = 15,
    # diff_alpha = 0.1,
    # ## delta_alpha = 0.5,
    # ## delta_alpha = 3.5,
    # delta_alpha = 3.75,
    diff_alpha = -0.28,
    delta_alpha = 1.4,
    # MAX_EFFECT = exp(-2.5),
    MAX_EFFECT = exp(-3),
    ## h = 4,
    h = 3,
    ## beta = exp(-2),
    beta = exp(-1.8),
    log_k = -1,
    nComp_sim = length(drug_regimen),
    T_nmblast = 5,
    T_retic = 5,
    T_RBC_max = 120,
    T_transit_steady_state = 3.5,
    G6PD_initial = 1,
    G6PD_decay_rate = 0.01,
    mu_death = log(0.31),
    sigma_death = 70
  )

  p <- ggplot(
    results |> pivot_longer(! time),
    aes(time, value)
  ) +
    geom_line() +
    geom_point(data = clin_data$all) +
    expand_limits(y = 0) +
    facet_wrap(~ name, scales = "free_y", ncol = 1)

  cat("Writing", plot_file, "...\n")
  ggsave(plot_file, p, width = 15, height = 25, unit = "cm")
}


inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

pow <- function(base, power) {
  base ** power
}

real2int <- function(x, min_val, max_val) {
  y <- round(x)
  for (n in min_val:max_val) {
    if (n >= y) {
      return(n)
    }
  }
  return(max_val)
}

dose_response <- function(dose, MAX_EFFECT, h, beta) {
  MAX_EFFECT * pow(dose, h) / (pow(dose, h) + pow(beta, h))
}

compute_rho <- function(Hb, Hb_delta, Hb_star, diff_alpha, delta_alpha) {
  exp(diff_alpha * (Hb - Hb_star)) * exp(delta_alpha * Hb_delta)
  # 1 + exp(diff_alpha * (Hb - Hb_star)) * exp(delta_alpha * Hb_delta)
  # max <- 5
  # max - (max - 1) * inv_logit(diff_alpha * (delta_alpha + Hb - Hb_star))
}

compute_transit_time <- function(Hb, Hb_star, T_transit_ss, log_k) {
  transit <- T_transit_ss * exp(- (Hb_star - Hb) * exp(log_k))
  min(transit, T_transit_ss)
}

CountRetics <- function(transit, reticulocytes) {
  Total_Retics <- 0
  T_retic <- length(reticulocytes)
  i <- T_retic
  transit_idx <- ceiling(transit)
  # NOTE: indexing starts at 1, not 0, and only goes up to 5.
  transit_int <- real2int(transit_idx, 1, 5)
  retic_rem <- ceiling(transit) - transit
  while(i > transit + 1) {
    Total_Retics <- reticulocytes[i]
    i <- i - 1
  }
  Total_Retics + retic_rem * reticulocytes[transit_int]
}

forwardsim <- function(drug_regimen, Hb_star, diff_alpha, delta_alpha,
                       MAX_EFFECT, h, beta, log_k, nComp_sim, T_nmblast,
                       T_retic, T_RBC_max, T_transit_steady_state,
                       G6PD_initial, G6PD_decay_rate, mu_death, sigma_death) {
  Hb_delta <- 0
  lambda <- 1000

  normoblasts <- array(lambda, dim = T_nmblast)
  reticulocytes <- array(lambda, dim = T_retic)
  erythrocytes <- array(lambda, dim = T_RBC_max)
  erythrocytes_G6PD <- array(lambda, dim = T_RBC_max)

  Hb <- array(0, dim = nComp_sim)
  rho <- array(0, dim = nComp_sim)
  retic_percent <- array(0, dim = nComp_sim)
  drug_effect <- array(0, dim = nComp_sim)
  max_eryth_age <- array(0, dim = nComp_sim)

  # Initial state (t = 1).
  t <- 1
  Hb[t] <- Hb_star
  transit <- compute_transit_time(Hb[t], Hb_star, T_transit_steady_state,
                                  log_k)
  rho[t] <- compute_rho(Hb[t], Hb_delta, Hb_star, diff_alpha, delta_alpha)
  drug_effect[t] <- 0
  erythrocytes_G6PD[1] <- G6PD_initial
  for (i in 2:T_RBC_max) {
    erythrocytes_G6PD[i] <- erythrocytes_G6PD[i - 1] *exp(-G6PD_decay_rate)
    erythrocytes[i] <- erythrocytes[i - 1] * inv_logit(
      (log(erythrocytes_G6PD[i - 1]) - mu_death) * sigma_death
    )
  }
  Total_Retics <- CountRetics(transit, reticulocytes)
  Total_Eryths <- sum(erythrocytes)
  C_0 <- Total_Retics + Total_Eryths
  retic_percent[t] <- 100 * Total_Retics / C_0
  max_eryth_age[t] <- sum(erythrocytes >= 0.5)

  # Simulate forward for t > 1.
  for (t in 2:nComp_sim) {
    if (t > 2) {
      Hb_delta <- Hb[t - 2] - Hb[t - 1]
    }
    rho[t] <- compute_rho(Hb[t - 1], Hb_delta, Hb_star,
                          diff_alpha, delta_alpha)
    transit <- compute_transit_time(Hb[t - 1], Hb_star, T_transit_steady_state,
                                    log_k)
    prev_nbs <- normoblasts
    normoblasts[1] <- rho[t] * lambda
    for (i in 2:T_nmblast) {
      normoblasts[i] <- prev_nbs[i - 1]
    }

    prev_retics <- reticulocytes
    reticulocytes[1] <- prev_nbs[T_nmblast]
    for (i in 2:T_retic) {
      reticulocytes[i] <- prev_retics[i - 1]
    }

    drug_effect[t] <- dose_response(drug_regimen[t - 1], MAX_EFFECT, h, beta)
    prev_eryths <- erythrocytes
    prev_G6PD <- erythrocytes_G6PD
    erythrocytes[1] <- prev_retics[T_retic]
    erythrocytes_G6PD[1] <- G6PD_initial
    for (i in 2:T_RBC_max) {
      erythrocytes_G6PD[i] <-prev_G6PD[i - 1] *
        exp(-G6PD_decay_rate - drug_effect[t])
      erythrocytes[i] <- prev_eryths[i - 1] *
        inv_logit((log(erythrocytes_G6PD[i - 1]) - mu_death) * sigma_death)
    }

    Total_Retics <- CountRetics(transit, reticulocytes)
    Total_Eryths <- sum(erythrocytes)
    C_t <- Total_Retics + Total_Eryths
    retic_percent[t] <- 100 * Total_Retics / C_t
    Hb[t] <- Hb_star * C_t / C_0
    max_eryth_age[t] <- sum(erythrocytes >= 0.5)

    if (t == 2) {
      cat("t = 2, Hb = ", Hb[2], "\n")
    }
  }

  data.frame(
    time = 1:nComp_sim,
    Hb = Hb,
    retic_percent = retic_percent,
    rho = rho,
    drug_effect = drug_effect,
    max_eryth_age  = max_eryth_age
  )
}


load_job_data <- function(job_i) {
  # NOTE: run `make_stan_fitting_datasets.R` to generate this file.
  stan_input_data_file <- file.path("Rout", "stan_data_list.RData")
  if (! file.exists(stan_input_data_file)) {
    stop("Missing file: ", stan_input_data_file)
  }

  data_env <- new.env()
  load(stan_input_data_file, envir = data_env)

  num_jobs <- length(data_env$dat_stan_list)
  if (job_i < 1 || job_i > num_jobs) {
    stop("Job number must be between 1 and ", num_jobs)
  }

  job_data <- data_env$dat_stan_list[[job_i]]
  job_data$log_MAX_EFFECT_prior_mean <- -2.5
  job_data$log_MAX_EFFECT_prior_sigma <- 0.5
  job_data$log_beta_mean <- -1.5
  job_data$log_k_sigma <- 0.25

  job_data
}


extract_drug_regimen <- function(job_data) {
  job_data$drug_regimen[
    job_data$ind_start_regimen[1]:job_data$ind_end_regimen[1]
  ]
}


extract_clinical_data <- function(job_data, experiment = 1) {
  ix <- experiment
  data_hb_cbc <- job_data$Hb_CBC[
    job_data$ind_start_Hb_CBC[ix]:job_data$ind_end_Hb_CBC[ix]
  ]
  time_hb_cbc <- job_data$t_sim_CBC_Hb[
    job_data$ind_start_Hb_CBC[ix]:job_data$ind_end_Hb_CBC[ix]
  ]
  data_hb_hc <-job_data$Hb_Haemocue[
    job_data$ind_start_Hb_hemocue[ix]:job_data$ind_end_Hb_hemocue[ix]
  ]
  time_hb_hc <-job_data$t_sim_hemocue_Hb[
    job_data$ind_start_Hb_hemocue[ix]:job_data$ind_end_Hb_hemocue[ix]
  ]
  data_retic <- job_data$Retic_data[
    job_data$ind_start_retic[ix]:job_data$ind_end_retic[ix]
  ]
  time_retic <- job_data$t_sim_retic[
    job_data$ind_start_retic[ix]:job_data$ind_end_retic[ix]
  ]

  df_retic <- data.frame(
    time = time_retic,
    value = data_retic
  )
  df_hb <- data.frame(
    time = c(time_hb_cbc, time_hb_hc),
    value = c(data_hb_cbc, data_hb_hc),
    kind = c(
      rep("CBC", length(data_hb_cbc)),
      rep("Haemocue", length(data_hb_hc))
    )
  )

  df_all_data <- rbind(
    df_retic |>
      mutate(name = "retic_percent"),
    df_hb[, c("time", "value")] |>
      mutate(name = "Hb")
  )

  list(retic = df_retic, hb = df_hb, all = df_all_data)
}


extract_all_data <- function() {
  all_exps <- load_job_data(1)
  dfs <- list()
  for (i in seq_len(all_exps$N_experiment)) {
    df <- extract_clinical_data(all_exps, i)$all |> mutate(experiment = i)
    dfs[[length(dfs) + 1]] <- df
  }
  df_out <- do.call(rbind, dfs)
}


# NOTE: run the script
main()
