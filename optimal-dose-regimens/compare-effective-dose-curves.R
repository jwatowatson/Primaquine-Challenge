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

  # NOTE: arguments are (drug_regimen, nComp_sim, weights, K_weights)
  # We only need to draw samples for the dose_weights array.
  compute_effective_dose <- settings$get_effective_dose_function()
  fit <- settings$load_fit()
  dose_weights_list <- get_dose_weights(fit)

  nComp_sim <- 29
  n_draws <- 1000
  draw_ixs <- seq(
    from = 1,
    to = length(dose_weights_list),
    length.out = n_draws
  ) |> as.integer()

  # Load the dose regimens administered in the ascending-dose study.
  df_asc <- load_ascending_dose_regimens(settings$weight_kg) |>
    select(ID, day, dosemgkg_adjusted) |>
    rename(dose_mgkg = dosemgkg_adjusted) |>
    mutate(optimal = FALSE)

  # Calculate effective doses for each administered regimen, and retain only
  # the daily lower and upper bounds for each regimen.
  df_eff <- get_effective_patient_doses(
    df_asc,
    dose_weights_list[draw_ixs],
    nComp_sim,
    compute_effective_dose
  ) |>
    group_by(ID, day) |>
    summarise(
      lower_mgkg = min(effective_mgkg),
      upper_mgkg = max(effective_mgkg),
      .groups = "drop"
    ) |>
    filter(
      upper_mgkg > 0
    ) |>
    mutate(ID_num = substring(ID, 6) |> as.integer())

  # Load the optimal dose regimens for a small number of thresholds.
  thresholds <- seq(0.5, 1.0, by = 0.1)
  df_optimal <- load_optimal_dose_regimens(settings, thresholds)

  # Calculate effective doses for each optimal regimen, and retain only the
  # daily lower and upper bounds for each regimen.
  df_opt <- get_effective_optimal_doses(
    df_optimal,
    dose_weights_list[draw_ixs],
    nComp_sim,
    compute_effective_dose
  ) |>
    group_by(duration, day) |>
    summarise(
      lower_mgkg = min(effective_mgkg),
      upper_mgkg = max(effective_mgkg),
      .groups = "drop"
    ) |>
    filter(
      upper_mgkg > 0
    )

  p_cmp_eff <- ggplot() +
    geom_ribbon(
      aes(day, ymin = lower_mgkg, ymax = upper_mgkg, fill = duration),
      df_opt,
    ) +
    geom_ribbon(
      aes(day, ymin = lower_mgkg, ymax = upper_mgkg),
      df_eff,
      fill = "#9f9f9f"
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    scale_fill_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    scale_x_continuous(
      "Day",
      breaks = 7 * 0:4
    ) +
    ylab(
      "Effective Dose (mg/kg)"
    ) +
    facet_wrap(
      ~ ID_num,
      ncol = 4,
      labeller = as_labeller(function(x) paste("ADPQ", x))
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.01),
      legend.justification = c(1, 0)
    )

  ggsave(
    "compare-effective-dose-curves.png",
    p_cmp_eff,
    width = 6,
    height = 9
  )
}


load_optimal_dose_regimens <- function(settings, inspect_thresholds) {
  optimal_regimens <- list()
  for (scenario in settings$scenarios) {
    df_results <- readRDS(scenario$results_file) |>
      mutate(threshold = - threshold) |>
      filter(threshold %in% inspect_thresholds)

    all_regimens <- read_csv(
      scenario$regimens_file,
      show_col_types = FALSE
    ) |>
      mutate(regimen_ix = row_number())

    df_optimal <- df_results |>
      left_join(all_regimens, by = "regimen_ix") |>
      mutate(duration = paste(scenario$duration, "days"))

    optimal_regimens[[length(optimal_regimens) + 1]] <- df_optimal
  }

  # Retain the (unique) optimal dose regimens.
  bind_rows(optimal_regimens) |>
    pivot_longer(
      starts_with("day_"),
      names_prefix = "day_",
      names_to = "day",
      names_transform = as.integer,
      values_to = "dose"
    ) |>
    replace_na(list(dose = 0)) |>
    select(duration, regimen, day, dose) |>
    unique() |>
    mutate(
      dose_mgkg = dose * settings$dose_unit_mg / settings$weight_kg,
      optimal = TRUE
    )
}


load_ascending_dose_regimens <- function(reference_weight_kg) {
  data_env <- new.env()
  load(file.path("..", "Data", "RBC_model_data.RData"), envir = data_env)
  data_env$PQdat |>
    filter(study == "Part1") |>
    mutate(day = Study_Day + 1) |>
    select(ID, day, dose, dosemgkg, weight) |>
    mutate(ID_num = substring(ID, 6) |> as.integer()) |>
    filter(dosemgkg > 0) |>
    # NOTE: calculate adjusted doses for the reference body weight.
    mutate(dosemgkg_adjusted = dosemgkg * weight / reference_weight_kg)
}


get_dose_weights <- function(fit) {
  dose_weights <- fit$draws("dose_weights", format = "df") |>
    as_tibble() |>
    select(starts_with("dose_weights["))

  dose_weights_mat <- dose_weights |>
    as.matrix() |>
    unname()

  purrr::map(
    seq_len(nrow(dose_weights_mat)),
    function(i) dose_weights_mat[i, ]
  )
}


get_effective_patient_doses <- function(
  df_asc, dose_weights_list, nComp_sim, fun
) {
  patient_IDs <- df_asc |> pull(ID) |> unique()
  patient_regimens <- lapply(
    patient_IDs,
    function(patient_ID) {
      regimen <- df_asc |> filter(ID == patient_ID) |> pull(dose_mgkg)
      c(regimen, rep(0, nComp_sim - length(regimen)))
    }
  )

  effective_tbls <- list()
  for (ix in seq_along(patient_IDs)) {
    ID <- patient_IDs[ix]
    regimen <- patient_regimens[[ix]]
    for (weights_ix in seq_along(dose_weights_list)) {
      weights <- dose_weights_list[[weights_ix]]
      effective_mgkg <- fun(
        regimen, nComp_sim, weights, length(weights)
      ) |> as.numeric()

      tbl <- tibble(
        day = 1:nComp_sim,
        effective_mgkg = effective_mgkg,
        ID = ID,
        sample = weights_ix
      )

      effective_tbls[[length(effective_tbls) + 1]] <- tbl
    }
  }

  bind_rows(effective_tbls)
}


get_effective_optimal_doses <- function(
  df_optimal, dose_weights_list, nComp_sim, fun
) {
  optimal_tbls <- list()
  for (durn in unique(df_optimal$duration)) {

    df_subset <- df_optimal |>
      filter(duration == durn)

    for (reg in unique(df_subset$regimen)) {

      regimen <- df_subset |>
        filter(regimen == reg) |>
        pull(dose_mgkg)

      regimen <- c(regimen, rep(0, nComp_sim - length(regimen)))

      for (weights_ix in seq_along(dose_weights_list)) {
        weights <- dose_weights_list[[weights_ix]]
        effective_mgkg <- fun(
          regimen, nComp_sim, weights, length(weights)
        ) |> as.numeric()

        tbl <- tibble(
          day = 1:nComp_sim,
          effective_mgkg = effective_mgkg,
          duration = durn,
          regimen = reg,
          sample = weights_ix
        )

        optimal_tbls[[length(optimal_tbls) + 1]] <- tbl
      }
    }
  }

  bind_rows(optimal_tbls)
}


main()
