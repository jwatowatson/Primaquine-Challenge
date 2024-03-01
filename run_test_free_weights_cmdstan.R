#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(fit_only = TRUE)

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

  job_numbers <- c(1, 2, 3, 4)
  for (job_number in job_numbers) {
    job_data <- utils$create_job_data(job_number, max_dose_delay)

    # Fit the model to data
    fit <- utils$fit_model(
      model, job_data, chains = chains, warmup = samples, samples = samples
    )

    # Save the results.
    results_file <- file.path(
      "Rout", paste0(prefix, "_job", job_number, ".rds")
    )
    cat("Saving job #", job_number, " to ", results_file, "\n", sep = "")
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

call_main("run_test_free_weights_cmdstan.R", main)
