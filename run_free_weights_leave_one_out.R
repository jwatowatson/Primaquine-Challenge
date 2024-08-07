#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages()

  # Load the Stan model
  model_file <- file.path(
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- cmdstan_model(stan_file = model_file)

  # Number of chains & samples
  chains <- 4
  samples <- 1000

  # Allow for delayed dose effects for up to 9 days.
  max_dose_delay <- 9
  prefix <- paste0(
    "pop_fit_mechanistic_G6PD_cmdstan_max_delay_",
    max_dose_delay
  )

  # Define initial values for some of the model parameters.
  init_list <- utils$initial_parameter_values(chains)

  jobs <- utils$create_job_data_ascending_dose_leave_one_out(
    max_dose_delay, quiet = TRUE
  )

  for (job_ix in seq_along(jobs)) {
    job_data <- jobs[[job_ix]]

    # Fit the model to data, using the initial haemoglobin measurement (plus
    # random effect) for the initial state, when predicting the left-out
    # individual's response.
    fit <- utils$fit_model(
      model, job_data, chains = chains, warmup = samples, samples = samples,
      init = init_list
    )

    # Save the results.
    results_file <- file.path(
      "Rout", paste0(prefix, "_leave_one_out_", job_ix, ".rds")
    )
    cat("Saving job #", job_ix, " to ", results_file, "\n", sep = "")
    fit$save_object(file = results_file)
  }

  invisible(0)
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

call_main("run_free_weights_leave_one_out.R", main)
