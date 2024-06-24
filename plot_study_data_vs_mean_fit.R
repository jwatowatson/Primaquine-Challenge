#!/usr/bin/env -S Rscript --vanilla
#
# Plot the study data versus the mean model fit.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  # Load the data collected in each part of the study.
  df_ascending <- load_study_data("Part1")
  df_single <- load_study_data("Part2")

  # Load the mean fit for the ascending dose study.
  fit_asc <- load_ascending_dose_mean_fit(utils)

  # Load the mean prediction for the single dose study, using the posterior
  # from fitting the model to the ascending dose study.
  ids <- levels(df_single$ID2)
  fit_single <- load_single_dose_mean_prediction(ids)

  # NOTE: colour individuals by their cumulative dose at day 10.
  colour_scale <- scale_colour_distiller(
    NULL,
    guide = "none",
    palette = "Spectral"
  )
  x_axis_scale <-  scale_x_continuous(
    "Days since start of primaquine",
    breaks = c(0, 10, 20, 28)
  )

  p_cum_dose <- ggplot() +
    geom_line(
      aes(
        Study_Day, cum_dosemgkg,
        colour = colour_by,
        group = ID2
      ),
      df_ascending
    ) +
    geom_vline(
      xintercept = 10,
      linetype = "dashed",
      colour = "#9f9f9f"
    ) +
    colour_scale +
    x_axis_scale +
    scale_y_continuous(
      "Total primaquine dose (mg/kg)",
      breaks = 2 * (0:4)
    ) +
    theme_bw()

  p_hb <- ggplot() +
    geom_line(
      aes(
        Study_Day, Mean_hb,
        colour = colour_by,
        group = ID2
      ),
      df_ascending |>
        filter(! is.na(Mean_hb))
    ) +
    geom_line(
      aes(Time, Mean),
      fit_asc$hb,
      linewidth = 1
    ) +
    colour_scale +
    x_axis_scale +
    scale_y_continuous(
      "Haemoglbin (g/dL)",
      breaks = 2 * (5:8)
    ) +
    theme_bw()

  p_retic <- ggplot() +
    geom_line(
      aes(
        Study_Day, Mean_retic,
        colour = colour_by,
        group = ID2
      ),
      df_ascending |>
        filter(! is.na(Mean_retic))
    ) +
    geom_line(
      aes(Time, Mean),
      fit_asc$retic,
      linewidth = 1
    ) +
    colour_scale +
    x_axis_scale +
    scale_y_continuous(
      "Reticulocytes (%)",
      limits = c(0, NA),
      breaks = c(5, 10, 15)
    ) +
    theme_bw()

  grob_ascending <- gridExtra::arrangeGrob(
    p_cum_dose + ggtitle("a"),
    p_hb + ggtitle("b"),
    p_retic + ggtitle("c"),
    nrow = 1
  )

  x_axis_scale_single <-  scale_x_continuous(
    "Days since start of primaquine",
    breaks = 4 * (0:3)
  )

  p_single_hb <- ggplot() +
    geom_line(
      aes(
        Study_Day, Mean_hb,
        colour = colour_by,
        group = ID
      ),
      df_single |>
        filter(! is.na(Mean_hb))
    ) +
    geom_line(
      aes(Study_Day, Mean),
      fit_single$hb,
      linewidth = 1,
      linetype = "dashed"
    ) +
    colour_scale +
    x_axis_scale_single +
    scale_y_continuous(
      "Haemoglbin (g/dL)",
      breaks = 2 * (5:8)
    ) +
    coord_cartesian(
      xlim = c(0, 14)
    ) +
    theme_bw()

  p_single_retic <- ggplot() +
    geom_line(
      aes(
        Study_Day, Mean_retic,
        colour = colour_by,
        group = ID
      ),
      df_single |>
        filter(! is.na(Mean_retic))
    ) +
    geom_line(
      aes(Study_Day, Mean),
      fit_single$retic,
      linewidth = 1,
      linetype = "dashed"
    ) +
    colour_scale +
    x_axis_scale_single +
    scale_y_continuous(
      "Reticulocytes (%)",
      limits = c(0, NA),
      breaks = 0:7
    ) +
    coord_cartesian(
      xlim = c(0, 14)
    ) +
    theme_bw()

  grob_single <- gridExtra::arrangeGrob(
    p_single_hb + ggtitle("d"),
    p_single_retic + ggtitle("e"),
    nrow = 1
  )

  grob_both <- gridExtra::arrangeGrob(
    grob_ascending,
    grob_single,
    nrow = 2
  )

  ggsave(
    "study-data-vs-model.png",
    grob_both,
    width = 9,
    height = 6
  )

}


load_study_data <- function(study_name = c("Part1", "Part2")) {
  data_env <- new.env()
  load(file.path("Data", "RBC_model_data.RData"), envir = data_env)

  study_name <- match.arg(study_name)

  df_study <- data_env$PQdat |>
    filter(study == study_name) |>
    group_by(ID2) |>
    mutate(
      cum_dosemgkg = cumsum(dosemgkg)
    ) |>
    ungroup() |>
    mutate(ID2 = factor_patients_by_number(ID2))

  if (study_name == "Part1") {
    colour_on_day <- 10
  } else {
    colour_on_day <- 7
  }

  df_study |>
    filter(Study_Day == colour_on_day) |>
    rename(colour_by = cum_dosemgkg) |>
    select(ID2, colour_by) |>
    inner_join(df_study, by = "ID2")
}


factor_patients_by_number <- function(patient_ids) {
  unique_ids <- unique(patient_ids)
  # Strip the "ADPQ " prefix and convert to integers.
  patient_numbers <- as.integer(substring(unique_ids, 5))
  # Sort the patients by number, rather than alphabetically.
  patient_order <- unique_ids[order(patient_numbers)]
  # Return an ordered factor.
  factor(patient_ids, levels = patient_order, ordered = TRUE)
}


load_ascending_dose_mean_fit <- function(utils) {
  fit <- readRDS(file.path(
    "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
  ))

  K_weights <- fit$metadata()$stan_variable_sizes$dose_weights
  max_dose_delay <- K_weights - 1
  job_number <- 3
  job_data <- utils$create_job_data(job_number, max_dose_delay, quiet = TRUE)

  df_fit_hb <- utils$get_fit_Hb_values(fit, job_data) |>
    as_tibble() |>
    # Convert days from 1..29 to 0..28.
    mutate(Time = Time - 1) |>
    group_by(Time) |>
    summarise(Mean = mean(Value))

  df_fit_retic <- utils$get_fit_retic_pcnt_values(fit, job_data) |>
    as_tibble() |>
    # Convert days from 1..29 to 0..28.
    mutate(Time = Time - 1) |>
    group_by(Time) |>
    summarise(Mean = mean(Value))

  list(hb = df_fit_hb, retic = df_fit_retic)
}


load_single_dose_mean_prediction <- function(ids) {
  # NOTE: we can take the mean fit for each individual and calculate the mean
  # of means, because the sample sizes are identical for each individual.
  pred_single <- new.env()
  sys.source("plot_free_weights_predict_single_dose.R", envir = pred_single)
  output_dir <- "Rout"
  results_regex <- "^.*_predict_single_dose_(.+)\\.rds$"
  results_files <- list.files(
    path = output_dir,
    pattern = results_regex,
    full.names = TRUE
  )
  df_pred <- pred_single$collect_predictions(results_files, ids)

  df_fit_hb <- df_pred |>
    filter(measure == "Hb") |>
    # Convert days from 1..29 to 0..28.
    mutate(Study_Day = Study_Day - 1) |>
    group_by(Study_Day) |>
    summarise(Mean = mean(Mean))

  df_fit_retic <- df_pred |>
    filter(measure == "retic_percent") |>
    # Convert days from 1..29 to 0..28.
    mutate(Study_Day = Study_Day - 1) |>
    group_by(Study_Day) |>
    summarise(Mean = mean(Mean))

  list(hb = df_fit_hb, retic = df_fit_retic)
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

call_main("plot_study_data_vs_mean_fit.R", main)
