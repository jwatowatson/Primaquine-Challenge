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

  # Calculate the mean drop in haemoglobin for each regimen, and the day on
  # which the nadir occurs.
  mean_hb_drops <- model_results |>
    filter(measure == "Haemoglobin (g/dL)") |>
    group_by(duration, regimen, draw) |>
    summarise(
      hb_drop = max(value) - min(value),
      nadir = day[which.min(value)],
      .groups = "drop"
    ) |>
    group_by(duration, regimen) |>
    summarise(
      hbdrop_mean = mean(hb_drop),
      hbdrop_lower95 = quantile(hb_drop, 0.025),
      hbdrop_upper95 = quantile(hb_drop, 0.975),
      nadir_mean = mean(nadir),
      nadir_lower95 = quantile(nadir, 0.025),
      nadir_upper95 = quantile(nadir, 0.975),
      .groups = "drop"
    ) |>
    pivot_longer(
      !c(duration, regimen),
      names_to = c("name", "measure"),
      names_sep="_"
    ) |>
    pivot_wider(
      names_from = measure,
      values_from = value
    )

  print(mean_hb_drops)
}


main()
