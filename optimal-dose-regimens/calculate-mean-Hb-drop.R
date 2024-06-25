#!/usr/bin/env -S Rscript --vanilla

main <- function() {
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  identify <- new.env()
  sys.source("identify-optimal-regimens.R", envir = identify)
  forward_sim <- settings$get_forward_sim_function()

  # Identify the matrix column that contains the results for our chosen
  # threshold of 1.0 g/dL.
  chosen_threshold <- 1.0
  df_best <- identify$load_optimal_regimens(settings, chosen_threshold)

  # For each regimen, simulate the model for 1000 random individuals.
  model_results <- identify$predict_outcomes(
    settings,
    df_best$long,
    forward_sim = forward_sim,
    intervals = FALSE
  )

  # Calculate the mean drop in haemoglobin for each regimen.
  mean_hb_drops <- model_results |>
    filter(measure == "Haemoglobin (g/dL)") |>
    group_by(duration, regimen, draw) |>
    summarise(
      hb_drop = max(value) - min(value),
      .groups = "drop"
    ) |>
    group_by(duration, regimen) |>
    summarise(
      hb_drop = mean(hb_drop),
      .groups = "drop"
    )

  print(mean_hb_drops)
}


main()
