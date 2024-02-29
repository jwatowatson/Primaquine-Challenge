#!/usr/bin/env -S Rscript --vanilla
#
# Plot fits for the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  # Load required packages
  suppressPackageStartupMessages(library(cmdstanr))
  suppressPackageStartupMessages(library(posterior, warn.conflicts = FALSE))
  suppressPackageStartupMessages(library(bayesplot, warn.conflicts = FALSE))
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  suppressPackageStartupMessages(library(ggplot2))

  color_scheme_set("brightblue")

  plot_results_for_each_job()
}

plot_results_for_each_job <- function() {
  # Identify the results file for each completed fit.
  results_files <- list.files(
    path = "Rout",
    pattern = "pop_fit_free_weights_cmdstan_job.*\\.rds",
    full.names = TRUE
  )

  for (results_file in results_files) {
    # Extract the job number from the filename.
    job_number <- as.numeric(sub(
      ".*pop_fit_free_weights_cmdstan_job(.*)\\.rds",
      "\\1",
      results_file
    ))

    # Construct and save the results plots.
    plot_fit_results(results_file, job_number)
  }
}

load_fit_job_data <- function(fit, job_number) {
  data_env <- new.env()
  load(file.path("Data", "RBC_model_data.RData"), envir = data_env)
  PQdat <- data_env$PQdat

  job_env <- new.env()
  sys.source("master_functions.R", envir = job_env)
  sys.source("run_test_free_weights_cmdstan.R", envir = job_env)
  K_weights <- fit$metadata()$stan_variable_sizes$dose_weights
  max_dose_delay <- K_weights - 1

  job_env$create_job_data(job_number, PQdat, max_dose_delay, quiet = TRUE)
}

plot_fit_results <- function(results_file, job_number) {
  fit <- readRDS(results_file)

  # Produce a density plot of the dose-weighting parameters.
  dose_weight_hists <- mcmc_areas(fit$draws("dose_weights"))

  Y_hat1 <- fit$draws("Y_hat")
  job_data <- load_fit_job_data(fit, job_number)

  # Construct plots.
  plot_list <- plot_fit_vs_data(fit, job_data, Y_hat1)

  # Provide a description for each job.
  job_names <- c(
    "Both studies, common individuals",
    "Both studies, different individuals",
    "Ascending-dose only",
    "Single-dose only"
  )

  # Save each plot.
  for (plot_name in names(plot_list)) {
    plot <- plot_list[[plot_name]]
    out_file <- paste0(
      "free-weights-model-plot-job-",
      job_number,
      "-",
      plot_name,
      ".png"
    )
    cat("Writing", out_file, "...")
    png(out_file, width = 8, height = 16, units = "in", res = 150)
    # NOTE: add the job name and description to the plot.
    print(plot +
            ggtitle(paste0("Job ", job_number, ": ", job_names[job_number])) +
            theme(plot.title = element_text(hjust = 0.5)))
    invisible(dev.off())
    cat("\n")
  }
}

plot_fit_vs_data <- function(fit ,job_data, Y_hat1) {
  num_exps <- job_data$N_experiment

  num_samples <- dim(Y_hat1)[1]
  num_chains <- dim(Y_hat1)[2]

  # Extract each variable into an array [samples, chains, time].
  final_dim <- dim(Y_hat1)[3]
  num_values <- final_dim / 3
  Y_Hb <- Y_hat1[, , seq(1, final_dim, by = 3)]
  Y_retic_percent <- Y_hat1[, , seq(2, final_dim, by = 3)]
  Y_effective_dose <- Y_hat1[, , seq(3, final_dim, by = 3)]

  dfs_Hb <- list()
  dfs_retic <- list()
  dfs_eff_dose <- list()
  dfs_dose <- list()
  dfs_hdata <- list() # NOTE: Hb data
  dfs_rdata <- list() # NOTE: retic data

  for (id in 1:num_exps) {
    ixs_regimen <- job_data$ind_start_regimen[id]:job_data$ind_end_regimen[id]
    num_values <- length(ixs_regimen)

    # Extract the Hb values.
    # TODO: there's a CBC_correction thing we need to account for?
    # fit$draws("CBC_correction")
    Hb <- array(Y_Hb[, , ixs_regimen], c(num_samples, num_chains, num_values))
    df_Hb <- as.data.frame(ftable(Hb)) |>
      mutate(
        Experiment = !!id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs_Hb[[length(dfs_Hb) + 1]] <- df_Hb

    # Extract the reticulocyte values.
    retic <- array(Y_retic_percent[, , ixs_regimen],
                   c(num_samples, num_chains, num_values))
    df_retic <- as.data.frame(ftable(retic)) |>
      mutate(
        Experiment = !!id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs_retic[[length(dfs_retic) + 1]] <- df_retic

    # Extract the effective dose values.
    eff_dose <- array(Y_effective_dose[, , ixs_regimen],
                      c(num_samples, num_chains, num_values))
    df_eff_dose <- as.data.frame(ftable(eff_dose)) |>
      mutate(
        Experiment = !!id,
        Sample = as.numeric(Var1),
        Chain = as.numeric(Var2),
        Time = as.numeric(Var3),
        Value = Freq
      ) |>
      select(Experiment, Sample, Chain, Time, Value)
    dfs_eff_dose[[length(dfs_eff_dose) + 1]] <- df_eff_dose

    # Extract the dosing interval.
    dosing_ind <- which(job_data$drug_regimen[ixs_regimen] > 0)
    dosing_span <- c(head(dosing_ind, 1), tail(dosing_ind, 1) + 1)
    df_dose <- data.frame(
      Experiment = id,
      Start = dosing_span[1],
      Until = dosing_span[2]
    )
    dfs_dose[[length(dfs_dose) + 1]] <- df_dose

    # Extract the experimental haemoglobin data.
    ixs_haemocue <- job_data$ind_start_Hb_hemocue[id]:job_data$ind_end_Hb_hemocue[id]
    t_haemocue <- job_data$t_sim_hemocue_Hb[ixs_haemocue]
    v_haemocue <- job_data$Hb_Haemocue[ixs_haemocue]

    ixs_hb_CBC <- job_data$ind_start_Hb_CBC[id]:job_data$ind_end_Hb_CBC[id]
    t_hb_CBC <- job_data$t_sim_CBC_Hb[ixs_hb_CBC]
    v_hb_CBC <- job_data$Hb_CBC[ixs_hb_CBC]

    df_hdata <- data.frame(
      Experiment = id,
      Time = c(t_haemocue, t_hb_CBC),
      Value = c(v_haemocue, v_hb_CBC),
      Measurement = c(
        rep("Haemocue", length(t_haemocue)),
        rep("CBC", length(t_hb_CBC))
      )
    )
    dfs_hdata[[length(dfs_hdata) + 1]] <- df_hdata

    # Extract the experimental reticulocyte data.
    ixs_retic_data <- job_data$ind_start_retic[id]:job_data$ind_end_retic[id]
    t_retic_data <- job_data$t_sim_retic[ixs_retic_data]
    v_retic_data <- job_data$Retic_data[ixs_retic_data]

    df_rdata <- data.frame(
      Experiment = id,
      Time = t_retic_data,
      Value = v_retic_data
    )
    dfs_rdata[[length(dfs_rdata) + 1]] <- df_rdata
  }

  df_Hb <- bind_rows(dfs_Hb)
  df_mean_Hb <- df_Hb |>
    group_by(Experiment, Time) |>
    summarise(
      Mean = mean(Value),
      Lower = quantile(Value, probs = 0.1),
      Upper = quantile(Value, probs = 0.9),
      .groups = "drop"
    )

  df_retic <- bind_rows(dfs_retic)
  df_mean_retic <- df_retic |>
    group_by(Experiment, Time) |>
    summarise(
      Mean = mean(Value),
      Lower = quantile(Value, probs = 0.1),
      Upper = quantile(Value, probs = 0.9),
      .groups = "drop"
    )

  df_eff_dose <- bind_rows(dfs_eff_dose)
  df_mean_eff_dose <- df_eff_dose |>
    group_by(Experiment, Time) |>
    summarise(
      Mean = mean(Value),
      Lower = quantile(Value, probs = 0.1),
      Upper = quantile(Value, probs = 0.9),
      .groups = "drop"
    )

  df_dose <- bind_rows(dfs_dose)
  df_hdata <- bind_rows(dfs_hdata)
  df_rdata <- bind_rows(dfs_rdata)

  plot_Hb <- ggplot() +
    geom_rect(
      aes(xmin = Start, xmax = Until, ymin = -Inf, ymax = Inf),
      df_dose,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Time, ymin = Lower, ymax = Upper),
      df_mean_Hb,
      fill = "#9f9f9f"
    ) +
    geom_line(aes(Time, Mean), df_mean_Hb) +
    geom_point(aes(Time, Value, colour = Measurement), df_hdata) +
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

  plot_retic <- ggplot() +
    geom_rect(
      aes(xmin = Start, xmax = Until, ymin = -Inf, ymax = Inf),
      df_dose,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Time, ymin = Lower, ymax = Upper),
      df_mean_retic,
      fill = "#9f9f9f"
    ) +
    geom_line(aes(Time, Mean), df_mean_retic) +
    geom_point(aes(Time, Value), df_rdata) +
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

  pars <- c(
    "Hb_star",
    "T_E_star",
    "alpha_diff1",
    "alpha_delta1",
    "alpha_diff2",
    "alpha_delta2",
    "logit_alpha",
    "beta",
    "h",
    "dose_weights", # NOTE: dose_weights[1] to dose_weights[15]
    "theta_rand"  # NOTE: theta_rand[1,1] to theta_rand[26,8]
  )
  pars_draws_df <- as_draws_df(fit$draws(pars))
  # NOTE: can convert to vanilla data frame with as.data.frame(), and then
  # consider using pivot_longer().
  # How to handle the non-scalar parameters?
  # Need to extract indices and convert into new column(s), with default
  # values of 1 for all other variables.

  test_df <- pars_draws_df |>
    as.data.frame() |>
    pivot_longer(!c(.chain, .iteration, .draw))

  # Plot the delay in dose effect.
  dose_weights_df <- test_df |>
    filter(startsWith(name, "dose_weights")) |>
    mutate(delay = as.integer(sub("dose_weights\\[(.+)\\]", "\\1", name)) - 1)
  plot_dose_weights <- ggplot() +
    geom_violin(aes(factor(delay), value), dose_weights_df) +
    xlab("Delay (days)") +
    ylab("Dose weight")

  # Reduction in RBC lifespan (%).
  dose_response <- function(effective_dose, logit_alpha, h, beta) {
    numer <- effective_dose^h
    denom <- effective_dose^h + beta^h
    100 * gtools::inv.logit(logit_alpha) * numer / denom
  }

  dose_response_df <- cross_join(
    # Effective Primaquine dose (mg/kg).
    data.frame(effective_dose = seq(0, 1, length.out = 200)),
    pars_draws_df |>
      as.data.frame() |>
      select(c(logit_alpha, beta, h, .draw))
  ) |>
    mutate(response = dose_response(effective_dose, logit_alpha, beta, h))
  dose_response_mean_df <- dose_response_df |>
    group_by(effective_dose) |>
    summarise(
      Mean = mean(response),
      Lower = quantile(response, probs = 0.1),
      Upper = quantile(response, probs = 0.9),
      .groups = "drop"
    )

  # Extract the random effects for each subject and plot the mean
  # dose-response curve for each subject.
  subject_responses <- list()
  for (subject_id in seq_len(job_data$N_subject)) {
    # logit_alpha + theta_rand[, id, 6]
    # beta + theta_rand[, id, 7]
    logit_alpha_effect_name <- paste0("theta_rand[", subject_id, ",6]")
    beta_effect_name <- paste0("theta_rand[", subject_id, ",7]")
    logit_alpha_effect <- pars_draws_df[[logit_alpha_effect_name]]
    beta_effect <- pars_draws_df[[beta_effect_name]]
    subject_response <- cross_join(
      # Effective Primaquine dose (mg/kg).
      data.frame(effective_dose = seq(0, 1, length.out = 200)),
      pars_draws_df |>
        as.data.frame() |>
        select(c(logit_alpha, beta, h, .draw)) |>
        mutate(
          logit_alpha = logit_alpha + logit_alpha_effect,
          beta = beta * exp(beta_effect)
        )
    ) |>
      mutate(response = dose_response(effective_dose, logit_alpha, beta, h)) |>
      group_by(effective_dose) |>
      summarise(Mean = mean(response)) |>
      mutate(subject = !!subject_id)
    subject_responses[[length(subject_responses) + 1]] <- subject_response
  }
  subject_response_df <- bind_rows(subject_responses)

  plot_dose_response <- ggplot() +
    geom_ribbon(
      aes(effective_dose, ymin = Lower, ymax = Upper),
      dose_response_mean_df,
      fill = "#afafaf"
    ) +
    geom_line(
      aes(effective_dose, Mean, group = subject),
      subject_response_df
    ) +
    geom_line(
      aes(effective_dose, Mean),
      dose_response_mean_df,
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

  # Plot a violin for each dose-response parameter.
  response_pars_df <- pars_draws_df |>
    as.data.frame() |>
    select(c(logit_alpha, beta, h, .draw)) |>
    pivot_longer(! .draw)
  plot_response_pars <- ggplot() +
    geom_violin(aes(name, value), response_pars_df) +
    xlab(NULL) +
    ylab(NULL)

  # Plot steady-state RBC lifespan (T_E_star) and haemoglobin (Hb_star)
  steady_state_df <- pars_draws_df |>
    as.data.frame() |>
    select(c(Hb_star, T_E_star, .draw)) |>
    pivot_longer(! .draw)
  # Enforce consistent bounds across all jobs.
  steady_state_bounds <- data.frame(
    name = rep(c("Hb_star", "T_E_star"), each = 2),
    value = c(13, 16, 50, 110)
  )
  plot_steady_state <- ggplot() +
    geom_violin(aes(name, value), steady_state_df) +
    geom_blank(aes(name, value), steady_state_bounds) +
    xlab(NULL) +
    ylab(NULL) +
    facet_wrap(~ name, nrow = 1, scale = "free")

  list(Hb = plot_Hb, retic = plot_retic, dose_weights = plot_dose_weights,
       dose_response = plot_dose_response, response_pars = plot_response_pars,
       steady_state = plot_steady_state)
}


call_main <- function(script_name, main) {
  if (interactive()) {
    return()
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=.+", args, value = TRUE)
  if (length(file_arg) != 1) {
    warning("Found ", length(file_arg), " '--file=' arguments")
    return()
  }

  file_path <- substring(file_arg, 8)
  file_name <- basename(file_path)
  if (file_name == script_name) {
    main_args <- commandArgs(trailingOnly = TRUE)
    status <- main(main_args)
    if (! is.numeric(status)) {
      status <- 0
    }
    quit(status = status)
  }
}

call_main("plot_test_free_weights_cmdstan.R", main)
