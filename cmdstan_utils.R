#
# Collect common functions for fitting Stan models with CmdStan and plotting
# the results.
#

# Load packages required to fit Stan models with CmdStan.
load_packages <- function(fit_only = FALSE) {
  suppressPackageStartupMessages(library(cmdstanr))
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)

  if (! fit_only) {
    suppressPackageStartupMessages(library(posterior, warn.conflicts = FALSE))
    suppressPackageStartupMessages(library(bayesplot, warn.conflicts = FALSE))
    suppressPackageStartupMessages(library(ggplot2))

    color_scheme_set("brightblue")
  }
}


create_job_data <- function(job_number, max_dose_delay, quiet = FALSE) {
  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  PQdat <- data_env$PQdat

  master <- new.env()
  sys.source("master_functions.R", envir = master)


  K_weights <- max_dose_delay + 1

  if (job_number == 1) {
    # Individuals recruited in both studies have the same parameters.
    job_data <- master$make_stan_dataset(
      my_data = PQdat,
      K_weights = K_weights
    )
  } else if (job_number == 2) {
    # Individuals recruited in both studies can have different parameters.
    job_data <- master$make_stan_dataset(
      my_data = PQdat,
      ID_subject = "ID2",
      K_weights = K_weights
    )
  } else if (job_number == 3) {
    # Only the ascending-dose study.
    job_data <- master$make_stan_dataset(
      my_data = PQdat |> filter(study == "Part1"),
      # Generate predictions for an ascending-dose patient.
      data_pred = PQdat |> filter(ID == "PQ45 13"),
      K_weights = K_weights
    )
  } else if (job_number == 4) {
    # Only the single-dose study.
    job_data <- master$make_stan_dataset(
      my_data = PQdat |> filter(study == "Part2"),
      # Generate predictions for an ascending-dose patient, assuming that the
      # model parameters do not differ between dosing regimens.
      data_pred = PQdat |> filter(study == "Part1", ID == "ADPQ 22"),
      K_weights = K_weights
    )
  } else {
    stop("Invalid job number: ", job_number)
  }

  # Define the prior for the delayed-dose weights.
  job_data$prior_weights <- rep(1, job_data$K_weights)

  # NOTE: remove data that contain NA values.
  # Typically: Y_true_haemocue, Y_true_HbCBC, Y_true_Retic.
  for (name in names(job_data)) {
    if (any(is.na(job_data[[name]]))) {
      if (! quiet) {
        cat("Removing job data '", name, "'; it contains NAs\n", sep = "")
      }
      job_data[[name]] <- NULL
    }
  }

  job_data
}


initial_parameter_values <- function(chains) {
  values <- list()
  for (i in seq_len(chains)) {
    name <- paste0("chain", i)
    values[[name]] <- list(
      Hb_star = rnorm(1, mean = 15, sd = 0.25),
      diff_alpha = 0.1,
      delta_alpha = 0.5,
      log_beta = -2,
      h = 4,
      log_k = -3,
      sigma_CBC = 0.5,
      sigma_haemocue = 0.5,
      sigma_retic = 0.5,
      CBC_correction = rnorm(1),
      log_MAX_EFFECT = -2.5,
      log_G6PD_decay_rate= -3,
      sigma_death = 10
    )
  }
  values
}


fit_model <- function(
  model, job_data, chains = 4, warmup = 1000, samples = 1000, refresh = 10,
  max_treedepth = 10, seed = 12345, init = NULL
) {
  model$sample(
    data = job_data,
    chains = chains,
    parallel_chains = parallel::detectCores(),
    seed = seed,
    init = init,
    max_treedepth = max_treedepth,
    iter_warmup = warmup,
    iter_sampling = samples,
    refresh = refresh
  )
}


get_fit_Hb_values <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  Y_hat <- fit$draws("Y_hat")
  num_samples <- dim(Y_hat)[1]
  num_chains <- dim(Y_hat)[2]
  final_dim <- dim(Y_hat)[3]
  num_values <- final_dim / 3

  Y_Hb <- Y_hat[, , seq(1, final_dim, by = 3)]
  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]
    num_values <- length(ixs_regimen)

    # Extract the Hb values.
    # NOTE: need to account for `fit$draws("CBC_correction")`
    # See the likelihood calculations in the Stan models.
    Hb <- array(Y_Hb[, , ixs_regimen], c(num_samples, num_chains, num_values))
    df_exp <- as.data.frame(ftable(Hb)) |>
      mutate(
        Experiment = .env$id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_retic_pcnt_values <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  Y_hat <- fit$draws("Y_hat")
  num_samples <- dim(Y_hat)[1]
  num_chains <- dim(Y_hat)[2]
  final_dim <- dim(Y_hat)[3]
  num_values <- final_dim / 3

  Y_retic_percent <- Y_hat[, , seq(2, final_dim, by = 3)]
  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]
    num_values <- length(ixs_regimen)

    retic <- array(Y_retic_percent[, , ixs_regimen],
                   c(num_samples, num_chains, num_values))
    df_exp <- as.data.frame(ftable(retic)) |>
      mutate(
        Experiment = .env$id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_effective_dose_values <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  Y_hat <- fit$draws("Y_hat")
  num_samples <- dim(Y_hat)[1]
  num_chains <- dim(Y_hat)[2]
  final_dim <- dim(Y_hat)[3]
  num_values <- final_dim / 3

  Y_effective_dose <- Y_hat[, , seq(3, final_dim, by = 3)]
  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]
    num_values <- length(ixs_regimen)

    eff_dose <- array(Y_effective_dose[, , ixs_regimen],
                      c(num_samples, num_chains, num_values))
    df_exp <- as.data.frame(ftable(eff_dose)) |>
      mutate(
        Experiment = .env$id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_dosing_intervals <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]

    dosing_ind <- which(job_data$drug_regimen[ixs_regimen] > 0)
    dosing_span <- c(head(dosing_ind, 1), tail(dosing_ind, 1) + 1)
    df_exp <- data.frame(
      Experiment = id,
      Start = dosing_span[1],
      Until = dosing_span[2]
    )
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_data_haemoglobin <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]

    ixs_haemocue <- job_data$ind_start_Hb_hemocue[id]:job_data$ind_end_Hb_hemocue[id]
    t_haemocue <- job_data$t_sim_hemocue_Hb[ixs_haemocue]
    v_haemocue <- job_data$Hb_Haemocue[ixs_haemocue]

    ixs_hb_CBC <- job_data$ind_start_Hb_CBC[id]:job_data$ind_end_Hb_CBC[id]
    t_hb_CBC <- job_data$t_sim_CBC_Hb[ixs_hb_CBC]
    v_hb_CBC <- job_data$Hb_CBC[ixs_hb_CBC]

    df_exp <- data.frame(
      Experiment = id,
      Time = c(t_haemocue, t_hb_CBC),
      Value = c(v_haemocue, v_hb_CBC),
      Measurement = c(
        rep("Haemocue", length(t_haemocue)),
        rep("CBC", length(t_hb_CBC))
      )
    )
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_data_reticulocytes <- function(fit, job_data) {
  num_exps <- job_data$N_experiment

  dfs <- list()

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]

    ixs_retic_data <- job_data$ind_start_retic[id]:job_data$ind_end_retic[id]
    t_retic_data <- job_data$t_sim_retic[ixs_retic_data]
    v_retic_data <- job_data$Retic_data[ixs_retic_data]

    df_exp <- data.frame(
      Experiment = id,
      Time = t_retic_data,
      Value = v_retic_data
    )
    dfs[[length(dfs) + 1]] <- df_exp
  }

  bind_rows(dfs)
}


get_fit_draws_wide <- function(fit, variables) {
  fit$draws(variables) |>
    as_draws_df() |>
    as.data.frame()
}


get_fit_draws_long <- function(fit, variables) {
  get_fit_draws_wide(fit, variables) |>
    pivot_longer(!c(.chain, .iteration, .draw))
}


get_fit_dose_delay_weights <- function(fit) {
  get_fit_draws_long(fit, "dose_weights") |>
    mutate(delay = as.integer(sub("dose_weights\\[(.+)\\]", "\\1", name)) - 1)
}


mean_and_interval <- function(df, lower = 0.1, upper = 0.9) {
  need_cols <- c("Experiment", "Time", "Value")
  if (! all(need_cols %in% names(df))) {
    stop("mean_and_interval(): one or more required columns are missing")
  }
  df |>
    group_by(Experiment, Time) |>
    summarise(
      Mean = mean(Value),
      Lower = quantile(Value, probs = 0.1),
      Upper = quantile(Value, probs = 0.9),
      .groups = "drop"
    )
}


plot_haemoglobin <- function(df_mean, df_dose_intervals, df_exp_data) {
  ggplot() +
    geom_rect(
      aes(xmin = Start, xmax = Until, ymin = -Inf, ymax = Inf),
      df_dose_intervals,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Time, ymin = Lower, ymax = Upper),
      df_mean,
      fill = "#9f9f9f"
    ) +
    geom_line(aes(Time, Mean), df_mean) +
    geom_point(aes(Time, Value, colour = Measurement), df_exp_data) +
    scale_colour_brewer(NULL, palette = "Dark2") +
    xlab("Time (days)") +
    ylab("Hb") +
    facet_wrap(
      ~ Experiment,
      scales = "free_x",
      ncol = 4
    ) +
    theme(
      legend.position = "top",
      # NOTE: hide facet labels.
      strip.background = element_blank(),
      strip.text = element_blank()
    )
}


plot_reticulocytes <- function(df_mean, df_dose_intervals, df_exp_data) {
  ggplot() +
    geom_rect(
      aes(xmin = Start, xmax = Until, ymin = -Inf, ymax = Inf),
      df_dose_intervals,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Time, ymin = Lower, ymax = Upper),
      df_mean,
      fill = "#9f9f9f"
    ) +
    geom_line(aes(Time, Mean), df_mean) +
    geom_point(aes(Time, Value), df_exp_data) +
    xlab("Time (days)") +
    ylab("Reticulocyte (%)") +
    facet_wrap(
      ~ Experiment,
      scales = "free_x",
      ncol = 4
    ) +
    theme(
      legend.position = "top",
      # NOTE: hide facet labels.
      strip.background = element_blank(),
      strip.text = element_blank()
    )
}


plot_dose_delay_weights <- function(df_weights) {
  ggplot() +
    geom_violin(
      aes(factor(delay), value),
      df_weights
    ) +
    xlab("Delay (days)") +
    ylab("Dose weight")
}


plot_subject_dose_responses <- function(df_mean, df_subjects) {
  ggplot() +
    geom_ribbon(
      aes(effective_dose, ymin = Lower, ymax = Upper),
      df_mean,
      fill = "#afafaf"
    ) +
    geom_line(
      aes(effective_dose, Mean, group = subject),
      df_subjects
    ) +
    geom_line(
      aes(effective_dose, Mean),
      df_mean,
      colour = "#ef0000",
      linewidth = 1
    ) +
    scale_x_continuous(
      "Primaquine (mg/kg)",
      breaks = (0:5) / 5
    ) +
    scale_y_continuous(
      "Reduction in RBC lifespan (%)",
      breaks = (0:6) * 10,
      limits = c(0, 60)
    )
}
