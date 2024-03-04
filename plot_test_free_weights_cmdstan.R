#!/usr/bin/env -S Rscript --vanilla
#
# Plot fits for the RBC model with independent dose-weighting parameters.
#
# USAGE:
#
#     ./plot_test_free_weights_cmdstan.R <MAX_DELAY>
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages()

  if (length(args) == 0) {
    # Default value
    max_delay <- 9
  } else if (length(args) == 1) {
    max_delay <- as.integer(args[1])
  } else {
    stop("Invalid arguments: ", args)
  }

  plot_results_for_each_job(utils, max_delay)
}

plot_results_for_each_job <- function(utils, max_delay) {
  re_match <- paste0(
    "pop_fit_free_weights_cmdstan_max_delay_",
    max_delay,
    "_job.*\\.rds"
  )
  re_job <- paste0(
    ".*pop_fit_free_weights_cmdstan_max_delay_",
    max_delay,
    "_job(.*)\\.rds"
  )

  # Identify the results file for each completed fit.
  results_files <- list.files(
    path = "Rout",
    pattern = re_match,
    full.names = TRUE
  )

  if (length(results_files) == 0) {
    warning("No output files for max_delay = ", max_delay)
  }

  for (results_file in results_files) {
    # Extract the job number from the filename.
    job_number <- as.numeric(sub(re_job, "\\1", results_file))

    # Define a filename prefix for plot files.
    plot_prefix <- paste0(
      "free-weights-model-max-delay-",
      max_delay,
      "-plot-job-"
    )

    # Construct and save the results plots.
    plot_fit_results(utils, results_file, job_number, plot_prefix)
  }
}

load_fit_job_data <- function(utils, fit, job_number) {
  K_weights <- fit$metadata()$stan_variable_sizes$dose_weights
  max_dose_delay <- K_weights - 1
  utils$create_job_data(job_number, max_dose_delay, quiet = TRUE)
}

plot_fit_results <- function(utils, results_file, job_number, plot_prefix) {
  fit <- readRDS(results_file)

  # Produce a density plot of the dose-weighting parameters.
  dose_weight_hists <- mcmc_areas(fit$draws("dose_weights"))

  # Construct plots.
  job_data <- load_fit_job_data(utils, fit, job_number)
  plot_list <- plot_fit_vs_data(utils, fit, job_data)

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
    out_file <- paste0(plot_prefix, job_number, "-", plot_name, ".png")
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

plot_fit_vs_data <- function(utils, fit, job_data) {
  df_Hb <- utils$get_fit_Hb_values(fit, job_data)
  df_mean_Hb <- utils$mean_and_interval(df_Hb)

  df_retic <- utils$get_fit_retic_pcnt_values(fit, job_data)
  df_mean_retic <- utils$mean_and_interval(df_retic)

  df_eff_dose <- utils$get_fit_effective_dose_values(fit, job_data)
  df_mean_eff_dose <- utils$mean_and_interval(df_eff_dose)

  df_dose_intervals <- utils$get_fit_dosing_intervals(fit, job_data)
  df_hb_data <- utils$get_fit_data_haemoglobin(git, job_data)
  df_retic_data <- utils$get_fit_data_reticulocytes(git, job_data)

  df_dose_delay_weights <- utils$get_fit_dose_delay_weights(fit)

  # Reduction in RBC lifespan (%).
  dose_response <- function(effective_dose, logit_alpha, h, beta) {
    numer <- effective_dose^h
    denom <- effective_dose^h + beta^h
    100 * gtools::inv.logit(logit_alpha) * numer / denom
  }

  dose_response_df <- cross_join(
    # Effective Primaquine dose (mg/kg).
    data.frame(effective_dose = seq(0, 1, length.out = 200)),
    utils$get_fit_draws_wide(fit, c("logit_alpha", "beta", "h"))
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
    alpha_effect_var <- paste0("theta_rand[", subject_id, ",6]")
    beta_effect_var <- paste0("theta_rand[", subject_id, ",7]")
    subject_response <- cross_join(
      # Effective Primaquine dose (mg/kg).
      data.frame(effective_dose = seq(0, 1, length.out = 200)),
      utils$get_fit_draws_wide(
        fit, c("logit_alpha", "beta", "h", alpha_effect_var, beta_effect_var)
      ) |>
        rename(
          alpha_effect = !!alpha_effect_var,
          beta_effect = !!beta_effect_var
        ) |>
        mutate(
          logit_alpha = logit_alpha + alpha_effect,
          beta = beta * exp(beta_effect)
        )
    ) |>
      mutate(response = dose_response(effective_dose, logit_alpha, beta, h)) |>
      group_by(effective_dose) |>
      summarise(Mean = mean(response)) |>
      mutate(subject = !!subject_id)
    subject_responses[[length(subject_responses) + 1]] <- subject_response
  }
  df_subject_responses <- bind_rows(subject_responses)

  # Plot model outputs vs experimental data.
  plot_Hb <- utils$plot_haemoglobin(df_mean_Hb, df_dose_intervals, df_hb_data)
  plot_retic <- utils$plot_reticulocytes(
    df_mean_retic, df_dose_intervals, df_retic_data
  )

  # Plot the dose delay weights and dose responses.
  plot_dose_weights <- utils$plot_dose_delay_weights(df_dose_delay_weights)
  plot_dose_response <- utils$plot_subject_dose_responses(
    dose_response_mean_df, df_subject_responses
  )

  # Plot a violin for each dose-response parameter.
  response_pars_df <- utils$get_fit_draws_long(
    fit, c("logit_alpha", "beta", "h")
  )
  plot_response_pars <- ggplot() +
    geom_violin(aes(name, value), response_pars_df) +
    xlab(NULL) +
    ylab(NULL)

  # Plot steady-state RBC lifespan (T_E_star) and haemoglobin (Hb_star)
  steady_state_df <- utils$get_fit_draws_long(
    fit, c("Hb_star", "T_E_star")
  )
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
