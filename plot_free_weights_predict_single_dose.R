#!/usr/bin/env -S Rscript --vanilla
#
# Plot the results when fitting the RBC model with independent dose-weighting
# to the ascending-dose study, and predicting the effect of the single-dose
# regimens.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  output_dir <- "Rout"
  results_regex <- "^.*_predict_single_dose_(.+)\\.rds$"
  results_files <- list.files(
    path = output_dir,
    pattern = results_regex,
    full.names = TRUE
  )
  left_out_ixs <- as.integer(sub(results_regex, "\\1",
                                 basename(results_files)))

  truth_dfs <- collect_ground_truth()
  ids <- unique(truth_dfs$hb$ID2)[left_out_ixs]

  df_pred <- collect_predictions(results_files, ids)

  p_retic <- plot_reticulocyte_percent(truth_dfs, df_pred)
  p_hb <- plot_haemoglobin(truth_dfs, df_pred)

  retic_file <- "predict-single-dose-retic-percent.png"
  cat("Writing", retic_file, "...")
  png(retic_file, width = 8, height = 8, units = "in", res = 150)
  print(p_retic)
  invisible(dev.off())
  cat("\n")

  hb_file <- "predict-single-dose-haemoglobin.png"
  cat("Writing", hb_file, "...")
  png(hb_file, width = 8, height = 8, units = "in", res = 150)
  print(p_hb)
  invisible(dev.off())
  cat("\n")

  invisible(0)
}


plot_reticulocyte_percent <- function(truth_dfs, df_pred) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  ggplot() +
    geom_rect(
      aes(xmin = Start_Day, xmax = Final_Day + 1, ymin = -Inf, ymax = Inf),
      truth_dfs$regimen,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Study_Day, ymin = Lower, ymax = Upper),
      df_pred |> filter(measure == "retic_percent"),
      fill = blues[2]
    ) +
    geom_line(
      aes(Study_Day, Median),
      df_pred |> filter(measure == "retic_percent"),
      colour = blues[3]
    ) +
    geom_point(
      aes(Study_Day, value, colour = name),
      truth_dfs$retic
    ) +
    scale_x_continuous(breaks = scales::breaks_width(7)) +
    scale_colour_brewer(NULL, palette = "Dark2") +
    xlab("Day") +
    ylab("Reticuloctye (%)") +
    expand_limits(y = 0) +
    facet_wrap(~ ID2, scale = "fixed", ncol = 4) +
    theme(legend.position = "top")
}


plot_haemoglobin <- function(truth_dfs, df_pred) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  ggplot() +
    geom_rect(
      aes(xmin = Start_Day, xmax = Final_Day + 1, ymin = -Inf, ymax = Inf),
      truth_dfs$regimen,
      fill = "#efcfef"
    ) +
    geom_ribbon(
      aes(Study_Day, ymin = Lower, ymax = Upper),
      df_pred |> filter(measure == "Hb"),
      fill = blues[2]
    ) +
    geom_line(
      aes(Study_Day, Median),
      df_pred |> filter(measure == "Hb"),
      colour = blues[3]
    ) +
    geom_point(
      aes(Study_Day, value, colour = name),
      truth_dfs$hb
    ) +
    scale_x_continuous(breaks = scales::breaks_width(7)) +
    scale_colour_brewer(NULL, palette = "Dark2") +
    xlab("Day") +
    ylab("Haemoglobin (g/dL)") +
    facet_wrap(~ ID2, scale = "fixed", ncol = 4) +
    theme(legend.position = "top")
}


collect_predictions <- function(results_files, ids) {
  predictions <- list()

  for (ix in seq_along(results_files)) {
    fit <- readRDS(results_files[ix])

    # Extract the model predictions for the left-out individual.
    draws <- as_draws_df(fit$draws("Y_pred")) |>
      as.data.frame() |>
      pivot_longer(! starts_with(".")) |>
      mutate(
        measure = case_when(
          startsWith(name, "Y_pred[1,") ~ "Hb",
          startsWith(name, "Y_pred[2,") ~ "retic_percent",
          startsWith(name, "Y_pred[3,") ~ "effective_dose"
        ),
        Study_Day = as.integer(sub("Y_pred\\[.,(\\d+)\\]", "\\1", name)),
        ID2 = ids[ix]
      ) |>
      select(! name)

    # Calculate the mean, median, and 5%-95% intervals.
    intervals <- draws |>
      group_by(ID2, Study_Day, measure) |>
      summarise(
        Mean = mean(value),
        Median = median(value),
        Lower = quantile(value, probs = 0.05),
        Upper = quantile(value, probs = 0.95),
        .groups = "drop"
      )

    predictions[[length(predictions) + 1]] <- intervals
  }

  bind_rows(predictions)
}


collect_ground_truth <- function() {
  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  study_data <- data_env$PQdat |> filter(study == "Part2")
  unique_ids <- unique(study_data$ID2)

  # Extract the data collected from each individual.
  # TODO: do we need to apply a correction to the reticulocyte data?
  df_true_retic <- study_data |>
    select(ID2, Study_Day, CBC_retic, Manual_retic) |>
    rename_with(~ gsub("_retic$", "", .x)) |>
    pivot_longer(! c(ID2, Study_Day)) |>
    filter(! is.na(value))

  df_true_hb <- study_data |>
    select(ID2, Study_Day, Haemocue_hb, CBC_hb) |>
    rename_with(~ gsub("_hb$", "", .x)) |>
    pivot_longer(! c(ID2, Study_Day)) |>
    filter(! is.na(value))

  # Extract the first and last day of treatment.
  df_regimen <- study_data |>
    filter(dosemgkg > 0) |>
    group_by(ID2) |>
    summarise(
      Start_Day = min(Study_Day),
      Final_Day = max(Study_Day),
    )

  list(retic = df_true_retic, hb = df_true_hb, regimen = df_regimen)
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

call_main("plot_free_weights_predict_single_dose.R", main)
