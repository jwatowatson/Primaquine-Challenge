#!/usr/bin/env -S Rscript --vanilla

## Load required packages
# packs <- c("dplyr", "tidyverse", "rstan",  "foreach", "parallel")
# invisible(lapply(packs, require, character.only = TRUE))

# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(tidyr)

color_scheme_set("brightblue")

source("master_functions.R")


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
  # job_data$log_MAX_EFFECT_prior_mean <- -2.5
  # job_data$log_MAX_EFFECT_prior_sigma <- 0.5
  job_data$log_MAX_EFFECT_prior_mean <- -2.75
  job_data$log_MAX_EFFECT_prior_sigma <- 1
  # job_data$log_beta_mean <- -1.5
  job_data$log_beta_mean <- -1.8
  job_data$log_k_sigma <- 0.25

  job_data
}


run_model <- function(model_file, chains, iters, thin, job_data) {
  model <- cmdstan_model(stan_file = model_file)
  job_data$Y_true_haemocue <- NULL
  job_data$Y_true_HbCBC <- NULL
  job_data$Y_true_Retic <- NULL
  out <- model$sample(
    data = job_data,
    # iter = iters,
    refresh = 10,
    chains = chains,
    seed = job_i,
    max_treedepth = 10
  )
  out
}


inspect_results <- function(stanfit, job_data, job_out_file) {
  shinystan::launch_shinystan(stanfit)
  sampler_params <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  summary(do.call(rbind, sampler_params), digits = 2)
  pairs(
    stanfit,
    pars = c("log_MAX_EFFECT", "log_G6PD_decay_rate"),
    las = 1
  )

  pairs(
    stanfit,
    pars = c(
      "diff_alpha", "delta_alpha",
      "log_G6PD_decay_rate","log_MAX_EFFECT",
      "sigma_death",
      "Hb_star", "h","log_beta","log_k"
    ),
    las = 1
  )

  rstan::traceplot(
    stanfit,
    pars = c(
      "diff_alpha", "delta_alpha",
      "log_G6PD_decay_rate","log_MAX_EFFECT",
      "sigma_death",
      "Hb_star", "h","log_beta","log_k",
      "sigma_CBC",
      "sigma_haemocue",
      # "sigma_retic",
      "CBC_correction"
    ),
    inc_warmup = FALSE
  )

  summary(
    stanfit,
    pars = c(
      "diff_alpha", "delta_alpha",
      "log_G6PD_decay_rate", "log_MAX_EFFECT",
      "sigma_death", # "mu_death",
      "Hb_star", "h","log_beta","log_k",
      "sigma_CBC",
      "sigma_haemocue",
      "sigma_retic",
      "CBC_correction"
    )
  )$summary

  thetas <- rstan::extract(stanfit, pars = "Y_hat")$Y_hat

  # NOTE: dimensions differ?!?
  # thetas has dims [1000, 4, 903]
  # job_data has length 798
  plot(job_data$t_sim_hemocue_Hb,
       colMeans(thetas[, 1, ]))
  plot(job_data$t_sim_hemocue_Hb,
       job_data$Hb_Haemocue-colMeans(thetas[,1,]))

  ind_retic <- get_ind_retic(job_data)
  plot(job_data$t_sim_retic,
       colMeans(thetas[, 2, ind_retic]))

  plot(colMeans(thetas[, 3, ind_retic]))
  plot(colMeans(thetas[, 4, ind_retic]))

  thetas <- rstan::extract(stanfit, pars = "Y_pred")$Y_pred

  y_hat <- rstan::extract(stanfit, pars = "Y_hat")$Y_hat
  y_pred <- rstan::extract(stanfit, pars = "Y_pred")$Y_pred

  par(mfrow = c(1, 1), las = 1)
  plot(
    colMeans(y_hat[, 1, ]), # ind_retic]),
    type = "l",
    col = "blue",
    ylab = "Hb"
  )
  points(
    colMeans(y_pred[, 1, ], na.rm = TRUE),
    col = "red"
  )
  points(
    job_data$t_sim_hemocue_Hb,
    job_data$Hb_Haemocue,
    col = "black",
    pch = 15
  )
  points(
    job_data$t_sim_CBC_Hb,
    job_data$Hb_CBC,
    col = "black",
    pch = 8
  )

  png(file="mechanistic-RBC-model-job-6-fit.png", width=800, height=800)
  par(mfrow = c(2, 2), las = 1)
  plot(
    colMeans(y_hat[, 1, ]), # ind_retic]),
    type = "l",
    col = "blue",
    ylab = "Hb"
  )
  points(
    colMeans(y_pred[, 1, ], na.rm = TRUE),
    col = "red"
  )
  points(
    job_data$t_sim_hemocue_Hb,
    job_data$Hb_Haemocue,
    col = "black",
    pch = 15
  )
  points(
    job_data$t_sim_CBC_Hb,
    job_data$Hb_CBC,
    col = "black",
    pch = 8
  )

  plot(
    colMeans(y_hat[, 2, ]), # ind_retic]),
    type = "l",
    col = "blue",
    ylab = "Reticulocytes (%)"
  )
  points(
    colMeans(y_pred[, 2, ], na.rm = TRUE),
    col = "red"
  )
  points(
    job_data$t_sim_retic,
    job_data$Retic_data,
    col = "black",
    pch = 15
  )

  plot(
    colMeans(y_hat[, 3, ind_retic]),
    type = "l",
    col = "blue",
    ylab = "Rho"
  )
  points(
    colMeans(y_pred[, 3, ], na.rm = TRUE),
    col = "red"
  )

  plot(
    colMeans(y_hat[, 4, ind_retic]),
    type = "l",
    col = "blue",
    ylab = "Drug effect"
  )
  points(
    colMeans(y_pred[, 4, ], na.rm = TRUE),
    col = "red"
  )
  dev.off()

  library(tidyr)

  df_out_hb <- as.data.frame(thetas[, 1, ]) |>
    mutate(sample = seq_len(dim(thetas)[1])) |>
    pivot_longer(! sample, names_to = "time", names_prefix = "V") |>
    mutate(time = as.numeric(time))
  df_exp_hb <- data.frame(
    time = job_data$t_sim_hemocue_Hb,
    value = job_data$Hb_Haemocue
  )

  # Plot experimental Hb (red points) vs model estimates (lines).
  ggplot() +
    geom_line(
      aes(time, value, group = sample),
      data = df_out_hb |> filter(time < 32, value < 20),
      alpha = 0.01
    ) +
    geom_point(
      aes(time, value),
      data = df_exp_hb |> filter(time < 32),
      colour = "red"
    )

  par(mfrow = c(2, 2), las = 1)

  # NOTE: this diverges from the data!
  # Likely caused by a small number of trajectories with massive
  # values, since we're using colMeans
  plot(
    # colMeans(thetas[, 1, ], na.rm = TRUE),
    apply(thetas[, 1, ], 2, function(x) median(x, na.rm = TRUE)),
    type = "l",
    lwd = 3,
    # ylim = c(8, 16),
    ylab = "Hb"
  )
  points(
    job_data$t_sim_hemocue_Hb,
    job_data$Hb_Haemocue,
    col = "red",
    pch = 15
  )
  plot(
    # colMeans(thetas[, 2 ,], na.rm = TRUE),
    apply(thetas[, 2, ], 2, function(x) median(x, na.rm = TRUE)),
    type = "l",
    lwd = 3,
    ylim = c(0, 15),
    ylab = "Retic (%)"
  )
  points(
    job_data$t_sim_retic,
    job_data$Retic_data,
    col = "red",
    pch = 15
  )

  plot(
    # colMeans(thetas[, 3, ], na.rm = TRUE),
    apply(thetas[, 3, ], 2, function(x) median(x, na.rm = TRUE)),
    type = "l",
    lwd = 3,
    ylab = "rho"
  )

  plot(
    # colMeans(thetas[, 4, ], na.rm = TRUE),
    apply(thetas[, 4, ], 2, function(x) median(x, na.rm = TRUE)),
    type = "l",
    lwd = 3,
    ylab = "Drug effect"
  )

  save(out, file = job_out_file)
  cat(sprintf("Finished job %s\n", job_i))
}


# Fitting details: model, chains, iterations.
model_file <- file.path("Stan_models", "RBC_model_mechanistic_G6PD_rm.stan")
chains <- 1
iters <- 1000
thin <- chains
job_i <- 6
job_out_file <- file.path("Rout", paste0("job", job_i, ".RData"))

job_data <- load_job_data(job_i)
out_file <- paste0("mech_rgm_obj_job_", job_i, "_debug.rds")
force_run <- FALSE
if (force_run || ! file.exists(out_file)) {
  out <- run_model(model_file, chains, iters, thin, job_data)
  out$save_output_files()
  out$save_object(file = out_file)
} else {
  out <- readRDS(out_file)
}

print(out$diagnostic_summary())

# inspect_results(out, job_data, job_out_file)

# NOTE: Y_pred[3,20] and Y_pred[3,21] contain 2 & 1 NA values, respectively.
## draws <- out$draws(format = "df")
## mcmc_hist(draws[, 2:13])

stanfit <- rstan::read_stan_csv(out$output_files())
