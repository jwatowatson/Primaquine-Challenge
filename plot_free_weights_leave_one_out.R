#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  output_dir <- "Rout"
  results_files <- list.files(
    path = output_dir,
    pattern = ".*_leave_one_out_.*\\.rds",
    full.names = TRUE
  )

  truth_dfs <- collect_ground_truth()
  patient_ids <- unique(truth_dfs$hb$ID2)
  df_pred <- collect_predictions(results_files, patient_ids)

  p_retic <- plot_reticulocyte_percent(truth_dfs, df_pred)
  p_hb <- plot_haemoglobin(truth_dfs, df_pred)

  retic_file <- "leave-one-out-retic-percent.png"
  cat("Writing", retic_file, "...")
  png(retic_file, width = 8, height = 12, units = "in", res = 150)
  print(p_retic)
  invisible(dev.off())
  cat("\n")

  hb_file <- "leave-one-out-haemoglobin.png"
  cat("Writing", hb_file, "...")
  png(hb_file, width = 8, height = 12, units = "in", res = 150)
  print(p_hb)
  invisible(dev.off())
  cat("\n")

  invisible(0)
}


plot_reticulocyte_percent <- function(truth_dfs, df_pred) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  ggplot() +
    geom_rect(
      aes(xmin = -Inf, xmax = Study_Day, ymin = -Inf, ymax = Inf),
      truth_dfs$final_doses,
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
      aes(Study_Day, CBC_retic),
      truth_dfs$retic
    ) +
    xlab("Day") +
    ylab("Reticuloctye (%)") +
    facet_wrap(~ ID2, scale = "free_x", ncol = 4)
}


plot_haemoglobin <- function(truth_dfs, df_pred) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  ggplot() +
    geom_rect(
      aes(xmin = -Inf, xmax = Study_Day, ymin = -Inf, ymax = Inf),
      truth_dfs$final_doses,
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
    scale_colour_brewer(NULL, palette = "Dark2") +
    xlab("Day") +
    ylab("Haemoglobin (g/dL)") +
    facet_wrap(~ ID2, scale = "free_x", ncol = 4) +
    theme(
      legend.position = c(1, 0), # "top")
      legend.justification = c(1, 0)
    )
}


collect_predictions <- function(results_files, unique_ids) {
  predictions <- list()

  for (ix in seq_along(results_files)) {
    results_file <- results_files[ix]
    left_out <- as.integer(
      sub("^.*_leave_one_out_(.+)\\.rds$", "\\1", results_file)
    )
    fit <- readRDS(results_file)

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
        ID2 = unique_ids[[left_out]]
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
  data_ascending <- data_env$PQdat |> filter(study == "Part1")
  unique_ids <- unique(data_ascending$ID2)

  # Extract the data collected from each individual.
  df_true_retic <- data_ascending |>
    select(ID2, Study_Day, CBC_retic) |>
    filter(! is.na(CBC_retic))

  df_true_hb <- data_ascending |>
    select(ID2, Study_Day, Haemocue_hb, CBC_hb) |>
    pivot_longer(! c(ID2, Study_Day)) |>
    filter(! is.na(value))

  df_final_doses <- data_ascending |>
    filter(dosemgkg > 0) |>
    group_by(ID2) |>
    summarise(Study_Day = max(Study_Day))

  list(retic = df_true_retic, hb = df_true_hb, final_doses = df_final_doses)
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

call_main("plot_free_weights_leave_one_out.R", main)
