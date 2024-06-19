#!/usr/bin/env -S Rscript --vanilla
#
# Plot the response to the effective dose, for the ascending-dose study fit.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  # Retrieve the dose_response() function from the Stan model.
  model_file <- file.path(
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- utils$compile_model_with_exposed_functions(model_file)
  model_fns <- model$functions

  # Load the model fit.
  # NOTE: Job #3 is the fit against the ascending-dose study.
  results_file <- file.path(
    "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
  )
  fit <- readRDS(results_file)

  # Reduction in RBC lifespan (%).
  df_mean_response <- mean_dose_response(fit)
  df_subject_responses <- subject_responses(fit)

  # Extract the random effects for each subject and plot the mean
  # dose-response curve for each subject.
  p_dose_response <- utils$plot_subject_dose_responses(
    df_mean_response, df_subject_responses
  ) +
    theme_bw()

  # Save the plot.
  plot_file <- "ascending-dose-fit-dose-response.png"
  cat("Writing", plot_file, "...")
  png(plot_file, width = 4.5, height = 4.5, units = "in", res = 150)
  print(p_dose_response)
  invisible(dev.off())
  cat("\n")

  invisible(0)
}

mean_dose_response <- function(fit) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)

  # Reduction in RBC lifespan (%).
  dose_response_df <- cross_join(
    # Effective Primaquine dose (mg/kg).
    data.frame(effective_dose = seq(0, 1, length.out = 201)),
    utils$get_fit_draws_wide(fit, c("logit_alpha", "beta", "h"))
  ) |>
    mutate(
      response = 100 * purrr:::pmap_dbl(
        list(effective_dose, logit_alpha, h, beta),
        model_fns$dose_response
      )
    )

  dose_response_df |>
    group_by(effective_dose) |>
    summarise(
      Mean = mean(response),
      Lower = quantile(response, probs = 0.025),
      Upper = quantile(response, probs = 0.975),
      .groups = "drop"
    )
}

subject_responses <- function(fit) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)

  num_subjects <- 23
  subject_responses <- list()

  for (subject_id in seq(23)) {

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
      mutate(
        response = 100 * purrr:::pmap_dbl(
          list(effective_dose, logit_alpha, h, beta),
          model_fns$dose_response
        )
      ) |>
      group_by(effective_dose) |>
      summarise(Mean = mean(response)) |>
      mutate(subject = !!subject_id)
    subject_responses[[length(subject_responses) + 1]] <- subject_response
  }

  bind_rows(subject_responses)
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

call_main("plot_dose_response.R", main)
