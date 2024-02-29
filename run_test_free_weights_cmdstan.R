#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  load_required_stuff()

  # Load the Stan model
  model_file <- file.path(
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- cmdstan_model(stan_file = model_file)

  # Number of chains & samples
  chains <- 4
  samples <- 1000

  # Allow for delayed dose effects for up to 14 days.
  max_dose_delay <- 14

  # Prepare the job data.
  load("Data/RBC_model_data.RData")

  job_numbers <- c(1, 2, 3, 4)
  for (job_number in job_numbers) {
    job_data <- create_job_data(job_number, PQdat, max_dose_delay)

    results_file <- file.path(
      "Rout", paste0("pop_fit_free_weights_cmdstan_job", job_number, ".rds")
    )

    # Fit the model to data
    fit <- model$sample(
      data = job_data,
      chains = chains,
      parallel_chains = parallel::detectCores(),
      seed = 12345,
      max_treedepth = 10,
      iter_warmup = samples,
      iter_sampling = samples,
      refresh = 10
    )

    # Save the results.
    cat("Saving job #", job_number, " to ", results_file, "\n", sep = "")
    fit$save_object(file = results_file)

    # mcmc_areas(fit$draws("dose_weights"))
  }

  invisible(0)
}

load_required_stuff <- function() {
  # Load required packages
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(dplyr)
  library(tidyr)

  color_scheme_set("brightblue")

  source("master_functions.R")
}

create_job_data <- function(job_number, PQdat, max_dose_delay, quiet = FALSE) {
  K_weights <- max_dose_delay + 1

  if (job_number == 1) {
    # Individuals recruited in both studies have the same parameters.
    job_data <- make_stan_dataset(
      my_data = PQdat,
      K_weights = K_weights
    )
  } else if (job_number == 2) {
    # Individuals recruited in both studies can have different parameters.
    job_data <- make_stan_dataset(
      my_data = PQdat,
      ID_subject = "ID2",
      K_weights = K_weights
    )
  } else if (job_number == 3) {
    # Only the ascending-dose study.
    job_data <- make_stan_dataset(
      my_data = PQdat |> filter(study == "Part1"),
      # Generate predictions for an ascending-dose patient.
      data_pred = PQdat |> filter(ID == "PQ45 13"),
      K_weights = K_weights
    )
  } else if (job_number == 4) {
    # Only the single-dose study.
    job_data <- make_stan_dataset(
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
