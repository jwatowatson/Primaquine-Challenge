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

  thresholds <- - seq(from = 0.1, to = 2, by = 0.1)

  results_file <- "evaluate-thresholds-results.rds"

  if (! file.exists(results_file)) {

    cat("Collecting results ...\n")

    plot_tables <- lapply(
      settings$scenarios,
      function(scenario) {
        evaluate_regimens(
          scenario$chunk_regex,
          scenario$num_chunks,
          thresholds,
          scenario$matrix_file,
          scenario$results_file
        )
      }
    )

    plot_tbl <- lapply(
      seq_along(plot_tables),
      function(ix) {
        plot_tables[[ix]]$plot |>
          mutate(duration = paste(settings$scenarios[[ix]]$duration, "days"))
      }
    ) |> bind_rows()

    # Save the summary data to disk.
    saveRDS(plot_tbl, results_file, compress = "xz")
  } else {
    cat("Loading results from", results_file, "...\n")
    plot_tbl <- readRDS(results_file)
  }

  p_cmp <- ggplot(plot_tbl) +
    geom_line(
      aes(threshold, 100 * fraction_below, colour = duration)
    ) +
    geom_point(
      aes(threshold, 100 * fraction_below, colour = duration)
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    xlab("Hb drop threshold (g/dL)") +
    ylab("Individuals exceeding threshold (%)") +
    expand_limits(
      x = 0,
      y = c(0, 100)
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.93, 0.93),
      legend.justification = c(1, 1)
    )

  ggsave(
    "individuals-exceeding-threshold.png",
    p_cmp,
    width = 6,
    height = 5
  )

  # NOTE: have to use factor(threshold) in order to adjust bar widths.
  p_num_optimal <- ggplot(
    plot_tbl |>
      filter(threshold > 0.3, threshold < 1.3)
  ) +
    geom_col(
      aes(factor(threshold), num_optimal, fill = duration),
      position = position_dodge2(),
      width = 0.6
    ) +
    scale_x_discrete(
      breaks = as.character(seq(0.4, 1.4, by = 0.2))
    ) +
    scale_fill_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    xlab("Hb drop threshold (g/dL)") +
    ylab("Number of optimal dose regimens") +
    theme(
      legend.position = "top"
    )

  ggsave(
    "individuals-exceeding-threshold-num-optimal.png",
    p_num_optimal,
    width = 6,
    height = 5
  )

  # Plot optimal dose regimens for a small number of thresholds.
  inspect_thresholds <- seq(0.5, 1.0, by = 0.1)
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

  df_optimal <- bind_rows(optimal_regimens) |>
    pivot_longer(
      starts_with("day_"),
      names_prefix = "day_",
      names_to = "day",
      names_transform = as.integer,
      values_to = "dose"
    ) |>
    filter(dose > 0)

  p_optimal <- ggplot(df_optimal) +
    geom_step(
      aes(day, dose * 2.5 / 60, colour = duration, group = regimen_ix)
    ) +
    scale_x_continuous(
      "Day",
      breaks = c(0, 5, 10, 15),
      limits = c(0, 15)
    ) +
    scale_y_continuous(
      "Daily dose (mg/kg)",
      limits = c(0, 0.8)
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    facet_wrap(
      ~ threshold,
      labeller = labeller(threshold = ~ paste("<=", .x, "g/dL"))
    ) +
    theme(
      legend.position = "top"
    )

  ggsave(
    "individuals-exceeding-threshold-optimal-regimens.png",
    p_optimal,
    width = 6,
    height = 5
  )

  # Plot the optimal regimens against those from the ascending-dose study.
  # The optimal regimens are not sensitive to the choice of threshold.
  data_env <- new.env()
  load(file.path("..", "Data", "RBC_model_data.RData"), envir = data_env)
  ascending_study <- data_env$PQdat |>
    filter(study == "Part1") |>
    mutate(day = Study_Day + 1) |>
    select(ID, day, dose, dosemgkg, weight) |>
    filter(dosemgkg > 0) |>
    # NOTE: calculate adjusted doses for the reference body weight.
    mutate(dosemgkg_adjusted = dosemgkg * weight / settings$weight_kg)

  p_vs_study <- ggplot() +
    geom_step(
      aes(
        day, dose * 2.5 / 60,
        group = interaction(regimen_ix, duration),
        colour = duration
      ),
      df_optimal
    ) +
    geom_step(
      aes(day, dosemgkg_adjusted, group = ID),
      ascending_study
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    scale_x_continuous(
      "Day",
      breaks = c(1, 7, 14, 21),
      minor_breaks = 0:21,
      limits = c(0, 21)
    ) +
    scale_y_continuous(
      "Dose (mg/kg)",
      limits = c(0, 0.8)
    ) +
    facet_wrap(~ ID, ncol = 4) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.01),
      legend.justification = c(1, 0)
    )

  ggsave(
    "individuals-exceeding-threshold-optimal-vs-ascending-study.png",
    p_vs_study,
    width = 6,
    height = 9
  )
}


fraction_below_threshold <- function(hb_drops, threshold) {
  rowSums(hb_drops < threshold) / ncol(hb_drops)
}


# Returns a matrix of dimensions N_REGIMENS x N_THRESHOLDS.
get_fractions_below_thresholds <- function(regex, num_chunks, thresholds) {
  chunk_files <- list.files(pattern = regex)
  if (length(chunk_files) != num_chunks) {
    stop("Expected ", num_chunks, " files but found ", length(chunk_files))
  }

  chunk_numbers <- as.integer(sub(regex, "\\1", chunk_files))

  chunked_results <- lapply(
    order(chunk_numbers),
    function(ix) {
      cat("Reading", chunk_files[ix], "...\n")
      chunk <- readRDS(chunk_files[ix])
      sapply(
        thresholds,
        function(threshold) {
          fraction_below_threshold(chunk, threshold)
        }
      )
    }
  )

  do.call(rbind, chunked_results)
}


# Returns a data frame of the optimal regimen(s) for each threshold, with
# columns 'regimen', 'threshold', 'fraction_below', 'num_optimal', and
# 'first_optimal_ix'.
fractions_data_frame <- function(fractions_mat, thresholds) {
  column_names <- paste0("threshold_", thresholds)
  colnames(fractions_mat) <- column_names

  as_tibble(fractions_mat) |>
    mutate(regimen = seq_len(nrow(fractions_mat))) |>
    pivot_longer(
      starts_with("threshold_"),
      names_prefix = "threshold_",
      names_to = "threshold",
      values_to = "fraction_below"
    ) |>
    mutate(
      threshold = as.numeric(threshold)
    ) |>
    group_by(threshold) |>
    mutate(regimen_ix = row_number()) |>
    filter(fraction_below == min(fraction_below)) |>
    # NOTE: record the # of optimal regimens for each threshold, and the index
    # of the first optimal regimen.
    mutate(
      num_optimal = n(),
      first_optimal_ix = dplyr::first(regimen_ix)
    ) |>
    ungroup()
}


evaluate_regimens <- function(regex, num_chunks, thresholds, matrix_file, results_file) {
  fracs_mat <- get_fractions_below_thresholds(regex, num_chunks, thresholds)
  saveRDS(fracs_mat, matrix_file, compress = "xz")

  results_tbl <- fractions_data_frame(fracs_mat, thresholds)
  saveRDS(results_tbl, results_file, compress = "xz")

  plot_tbl <- results_tbl |>
    select(threshold, fraction_below, num_optimal, first_optimal_ix) |>
    unique() |>
    mutate(threshold = - threshold)

  list(results = results_tbl, plot = plot_tbl)
}


main()
